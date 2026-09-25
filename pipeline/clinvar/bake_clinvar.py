"""Bake exact, reference-checked ClinVar SNVs; never infer a classification.

The source is the current VCV XML, not the legacy ClinicalSignificance schema.
All searched IDs must arrive before replacing an asset. Raw responses are cached
with retrieval dates and checksums for audit/replay, separate from the app.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import http.client
import json
import os
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))
from pipeline.check_assets import positions, bases, translate  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import TARGETS  # noqa: E402

COMPLEMENT = str.maketrans('ACGT', 'TGCA')
SCOPE = 'reference-matched single nucleotide variants in the drawn gene'

# VCV records per efetch. NCBI drops a long run's connection now and then, and
# Biopython retries only the request, never the read of its body, so a batch is
# fetched, read and parsed as one attempt and retried whole, with backoff.
BATCH = 100
ATTEMPTS = 5
BACKOFF_SECONDS = 10


def project_records(root: ET.Element, target, mock: dict, impact: dict) -> tuple[list, dict]:
    """Match GRCh38 alleles to the exact drawn bases (including reverse genes).

    Read only direct SimpleAllele locations: Gene/Location spans whole genes,
    and assertions/haplotypes may contain other variants with other meanings.
    """
    local_by_genomic = {}
    for run in impact['runs']:
        for offset in range(run['length']):
            genomic = run['genomic'] + run['step'] * offset
            if genomic in local_by_genomic:
                raise ValueError('Ambiguous genomic mapping')
            local_by_genomic[genomic] = run['local'] + offset
    # The chromosome span the drawn gene covers. Inside it, the only bases the
    # runs leave out are the middles of introns drawn shortened (DMD, APP,
    # CFTR): in the gene, but not on the page, which is a different fact from
    # lying outside it.
    first, last = min(local_by_genomic), max(local_by_genomic)
    cds = positions(mock, mock['protein']['segments'])
    cds_index = {p: i for i, p in enumerate(cds)}
    dna = bases(mock, cds)
    if translate(dna).rstrip('*') != mock['protein']['translation']:
        raise ValueError('CDS differs from the protein')
    variants, excluded, seen = [], Counter(), set()
    for record in root.findall('VariationArchive'):
        variation_id = record.attrib['VariationID']
        if variation_id in seen:
            raise ValueError(f'Duplicate Variation ID {variation_id}')
        seen.add(variation_id)
        if record.findtext('RecordStatus') != 'current':
            excluded['not_current'] += 1
            continue
        allele = record.find('ClassifiedRecord/SimpleAllele')
        if allele is None or not any(g.get('Symbol') == target.gene for g in allele.findall('GeneList/Gene')):
            excluded['not_a_simple_allele_of_gene'] += 1
            continue
        locations = [loc for loc in allele.findall('Location/SequenceLocation')
                     if loc.get('Assembly') == 'GRCh38' and 'chr' + loc.get('Chr', '') == impact['chromosome']]
        if len(locations) != 1:
            excluded['no_unique_GRCh38_location'] += 1
            continue
        loc = locations[0]
        ref, alt = loc.get('referenceAlleleVCF', ''), loc.get('alternateAlleleVCF', '')
        if (allele.findtext('VariantType') != 'single nucleotide variant' or
                len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt or
                loc.get('start') != loc.get('stop') or loc.get('start') != loc.get('positionVCF')):
            excluded['not_a_single_base_substitution'] += 1
            continue
        genomic = int(loc.attrib['positionVCF'])
        local = local_by_genomic.get(genomic)
        if local is None:
            excluded['intron_not_drawn' if first <= genomic <= last else 'outside_drawn_gene'] += 1
            continue
        genomic_ref, genomic_alt = ref, alt
        if impact['complemented']:
            ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
        if bases(mock, [local]) != ref:
            excluded['reference_mismatch'] += 1
            continue
        classification = record.find('ClassifiedRecord/Classifications/GermlineClassification')
        if classification is None or not classification.findtext('Description'):
            excluded['no_germline_classification'] += 1
            continue
        conditions = []
        for rcv in record.findall('ClassifiedRecord/RCVList/RCVAccession'):
            germline = rcv.find('RCVClassifications/GermlineClassification')
            if germline is None:
                continue
            description = germline.find('Description')
            conditions.append({
                'accession': rcv.attrib['Accession'] + '.' + rcv.attrib['Version'],
                'names': [c.text for c in rcv.findall('ClassifiedConditionList/ClassifiedCondition') if c.text],
                'classification': germline.findtext('Description') or 'Not provided',
                'review_status': germline.findtext('ReviewStatus') or 'Not provided',
                'last_evaluated': description.get('DateLastEvaluated') if description is not None else None,
            })
        offset = cds_index.get(local)
        residue, protein_change, coding_change, consequence = None, None, None, None
        if offset is not None:
            coding_change = f'c.{offset + 1}{ref}>{alt}'
            codon_start = offset // 3 * 3
            triplet = dna[codon_start:codon_start + 3]
            aa = translate(triplet)
            changed = translate(triplet[:offset % 3] + alt + triplet[offset % 3 + 1:])
            if aa != '*':
                residue = offset // 3 + 1
                protein_change = f'p.{aa}{residue}{changed if changed != aa else "="}'
                consequence = ('synonymous' if changed == aa else 'start codon' if residue == 1
                               else 'stop gained' if changed == '*' else 'missense')
            else:
                consequence = 'stop retained' if changed == '*' else 'stop lost'
        variants.append({
            'variation_id': variation_id,
            'accession': record.attrib['Accession'] + '.' + record.attrib['Version'],
            'name': record.attrib['VariationName'],
            'position': local, 'genomic': genomic,
            'ref': ref, 'alt': alt, 'genomic_ref': genomic_ref, 'genomic_alt': genomic_alt,
            'residue': residue, 'protein_change': protein_change,
            'coding_change': coding_change, 'consequence': consequence,
            'classification': classification.findtext('Description'),
            'review_status': classification.findtext('ReviewStatus') or 'Not provided',
            'last_evaluated': classification.get('DateLastEvaluated'),
            'last_updated': record.get('DateLastUpdated'),
            'conditions': conditions,
            'collection_methods': sorted({e.text for e in record.findall('ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion/ObservedInList/ObservedIn/Method/MethodType') if e.text}),
        })
    return sorted(variants, key=lambda v: (v['position'], v['alt'], int(v['variation_id']))), dict(sorted(excluded.items()))


PLACEHOLDERS = {'not provided', 'not specified'}


def condition_traits(root: ET.Element, variants: list) -> dict:
    """Name each mapped record's conditions by the identifiers ClinVar gives them.

    An RCV names its condition and, where ClinVar has mapped it, its MedGen
    concept; the record's own ConditionList carries that concept's OMIM, MONDO
    and symbol. The symbol is the one OMIM cross-references for the condition's
    own entry, else ClinVar's preferred one. Definitions are deliberately not
    taken: MedGen attaches the WFS1 GeneReviews summary to type 2 diabetes.
    """
    mapped = {v['variation_id'] for v in variants}
    traits = {}
    for record in root.findall('VariationArchive'):
        if record.attrib['VariationID'] not in mapped:
            continue
        by_concept = {}
        for trait in record.findall('ClassifiedRecord/Classifications/GermlineClassification/ConditionList/TraitSet/Trait'):
            for xref in trait.findall('XRef'):
                if xref.get('DB') == 'MedGen':
                    by_concept[xref.get('ID')] = trait
        for rcv in record.findall('ClassifiedRecord/RCVList/RCVAccession'):
            if rcv.find('RCVClassifications/GermlineClassification') is None:
                continue
            for condition in rcv.findall('ClassifiedConditionList/ClassifiedCondition'):
                name = condition.text
                if not name or name in PLACEHOLDERS or condition.get('DB') != 'MedGen':
                    continue
                entry = {'medgen': condition.get('ID')}
                trait = by_concept.get(entry['medgen'])
                if trait is not None:
                    omim = [x.get('ID') for x in trait.findall('XRef') if x.get('DB') == 'OMIM']
                    own = next((i for i in omim if not i.startswith('PS')), omim[0] if omim else None)
                    mondo = [x.get('ID') for x in trait.findall('XRef') if x.get('DB') == 'MONDO']
                    # Records list a trait's symbols in no fixed order (type 2
                    # diabetes: T2D and NIDDM, both on OMIM 125853), so ties are
                    # broken by rank, never by position: preferred first, then
                    # the one fewest other sources claim, then the shorter.
                    symbols = sorted(
                        ((s.find('ElementValue').get('Type') != 'Preferred', len(s.findall('XRef')),
                          len(s.findtext('ElementValue')), s.findtext('ElementValue'),
                          {x.get('ID') for x in s.findall('XRef') if x.get('DB') == 'OMIM'})
                         for s in trait.findall('Symbol') if s.findtext('ElementValue')),
                        key=lambda s: s[:4])
                    symbol = next((s[3] for s in symbols if own is not None and own in s[4]),
                                  next((s[3] for s in symbols if not s[0]), None))
                    if symbol:
                        entry['symbol'] = symbol
                    if own:
                        entry['omim'] = own
                    if mondo:
                        entry['mondo'] = mondo[0]
                known = traits.setdefault(name, {})
                for key, value in entry.items():
                    if known.setdefault(key, value) != value:
                        raise ValueError(f'Condition {name!r} has conflicting identifiers')
    return dict(sorted(traits.items()))


def fetch_batch(efetch, ids: list[str], attempts: int = ATTEMPTS, pause=time.sleep) -> bytes:
    """One batch of VCV XML holding exactly the records asked for.

    A dropped connection, a truncated body and a response missing some of the
    records all look alike from here, and all get the same answer: the whole
    batch again, after a wait that doubles. The last failure is raised, which
    stops the bake before anything is written.
    """
    for attempt in range(1, attempts + 1):
        try:
            with efetch(db='clinvar', rettype='vcv', is_variationid='true', from_esearch='true',
                        id=','.join(ids)) as handle:
                raw = handle.read()
            if isinstance(raw, str):
                raw = raw.encode()
            parsed = ET.fromstring(raw)
            if parsed.tag != 'ClinVarResult-Set' or parsed.find('.//ERROR') is not None:
                raise ValueError('Unexpected ClinVar response')
            received = sorted(r.attrib['VariationID'] for r in parsed.findall('VariationArchive'))
            if received != sorted(ids):
                raise ValueError(f'Batch returned {len(received)} of the {len(ids)} records asked for')
            return raw
        except (OSError, http.client.HTTPException, ET.ParseError, KeyError, ValueError) as error:
            if attempt == attempts:
                raise
            wait = BACKOFF_SECONDS * 2 ** (attempt - 1)
            print(f'    batch of {len(ids)}: {error!r}; retrying in {wait}s', file=sys.stderr, flush=True)
            pause(wait)
    raise AssertionError('unreachable')


def fetch(gene: str, cache: Path, replay: bool) -> tuple[ET.Element, dict]:
    if replay:
        metadata = json.loads((cache / 'source.json').read_text())
        if metadata['query'] != f'{gene}[gene]':
            raise ValueError('Cached query differs')
    else:
        from Bio import Entrez
        Entrez.email = os.environ.get('NCBI_EMAIL', '')
        if not Entrez.email:
            raise ValueError('Set NCBI_EMAIL for NCBI Entrez requests')
        Entrez.tool = 'HelixPeek_ClinVar_baker'
        Entrez.api_key = os.environ.get('NCBI_API_KEY')
        ids, count = [], None
        for start in range(0, 1000000, 1000):
            with Entrez.esearch(db='clinvar', term=f'{gene}[gene]', retstart=start, retmax=1000, retmode='xml') as handle:
                result = Entrez.read(handle)
            current = int(result['Count'])
            if count is not None and current != count:
                raise ValueError('Search changed during pagination; retry')
            count = current
            ids.extend(map(str, result['IdList']))
            if len(ids) >= count:
                break
        if len(ids) != count or len(set(ids)) != count:
            raise ValueError('Incomplete/duplicate search IDs')
        cache.mkdir(parents=True, exist_ok=True)
        batches = []
        for start in range(0, len(ids), BATCH):
            raw = fetch_batch(Entrez.efetch, ids[start:start + BATCH])
            name = f'batch-{start // BATCH}.xml'
            (cache / name).write_bytes(raw)
            batches.append({'file': name, 'sha256': hashlib.sha256(raw).hexdigest()})
            print(f'    {min(start + BATCH, len(ids)):,} of {len(ids):,} records', file=sys.stderr, flush=True)
        metadata = {'query': f'{gene}[gene]', 'retrieved_at': datetime.now(timezone.utc).isoformat(),
                    'record_count': count, 'ids': ids, 'batches': batches}
        (cache / 'source.json').write_text(json.dumps(metadata, indent=2) + '\n')
    root = ET.Element('ClinVarResult-Set')
    for batch in metadata['batches']:
        raw = (cache / batch['file']).read_bytes()
        if hashlib.sha256(raw).hexdigest() != batch['sha256']:
            raise ValueError('Cached XML checksum mismatch')
        parsed = ET.fromstring(raw)
        if parsed.tag != 'ClinVarResult-Set' or parsed.find('.//ERROR') is not None:
            raise ValueError('Unexpected ClinVar response')
        root.extend(parsed.findall('VariationArchive'))
    received = [r.attrib['VariationID'] for r in root]
    if len(received) != metadata['record_count'] or set(received) != set(metadata['ids']):
        raise ValueError('Fetched records do not match all search IDs; asset not replaced')
    return root, metadata


def bake(target, replay: bool) -> None:
    mock = json.loads((DATA / target.mock_asset).read_text())
    impact = json.loads((DATA / target.impact_asset).read_text())
    if impact['gene'] != target.gene or impact['assembly'] != 'GRCh38' or impact['sequence'] != bases(mock, list(range(mock['location']['start'], mock['location']['end'] + 1))):
        raise ValueError('Impact mapping and gene record differ')
    root, source = fetch(target.gene, Path(__file__).parent / 'cache' / target.gene, replay)
    variants, excluded = project_records(root, target, mock, impact)
    if len(variants) + sum(excluded.values()) != source['record_count']:
        raise ValueError('Unaccounted records')
    asset = {'schema_version': 1, 'source': 'NCBI ClinVar', 'gene': target.gene,
             'accession': target.source.accession, 'assembly': 'GRCh38', 'chromosome': impact['chromosome'],
             'scope': SCOPE, 'start': impact['start'], 'sequence': impact['sequence'],
             'protein_sequence': mock['protein']['translation'], 'complemented': impact['complemented'],
             'runs': impact['runs'], 'retrieved_at': source['retrieved_at'], 'query': source['query'],
             'searched_records': source['record_count'], 'excluded': excluded,
             'source_batches': source['batches'], 'traits': condition_traits(root, variants),
             'variants': variants}
    path = DATA / f'assets/clinvar/{target.slug}_clinvar.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(asset, indent=2) + '\n')
    temporary.replace(path)
    print(f'{target.gene}: {len(variants)} mapped SNVs; {excluded}; {path.relative_to(DATA)}', flush=True)


def chosen(args) -> list:
    """The rows a run bakes: one by slug, or every row the table says has a
    snapshot — the same meaning `--all` has for the other bakers."""
    if args.all:
        return [t for t in TARGETS if t.clinvar_available]
    return [next(t for t in TARGETS if t.slug == args.target)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all', action='store_true', help='every row with clinvar_available')
    group.add_argument('--target', choices=[t.slug for t in TARGETS], help='one row, by slug')
    parser.add_argument('--replay', action='store_true', help='Use checksummed local XML; no network')
    args = parser.parse_args()
    failures = []
    for target in chosen(args):
        try:
            bake(target, args.replay)
        except ValueError as error:
            print(f'{target.gene}: FAILED: {error}', file=sys.stderr, flush=True)
            failures.append(target.slug)
    if failures:
        print(f'{len(failures)} failed: {", ".join(failures)}; their assets were not replaced', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
