"""Synthetic XML tests of the mapping boundary, independent of network access."""
import copy
import hashlib
import http.client
import json
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from types import SimpleNamespace

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))
from pipeline.clinvar.bake_clinvar import (  # noqa: E402
    chosen, condition_traits, fetch, fetch_batch, project_records,
)


def record(id='1', alt='C', classification='Conflicting classifications of pathogenicity', at=102):
    return ET.fromstring(f'''<VariationArchive VariationID="{id}" Accession="VCV00000000{id}" Version="2" VariationName="synthetic allele">
    <RecordStatus>current</RecordStatus><ClassifiedRecord><SimpleAllele>
      <GeneList><Gene Symbol="INS"><Location><SequenceLocation Assembly="GRCh38" Chr="11" start="1" stop="9"/></Location></Gene></GeneList>
      <VariantType>single nucleotide variant</VariantType><Location>
        <SequenceLocation Assembly="GRCh38" Chr="11" start="{at}" stop="{at}" positionVCF="{at}" referenceAlleleVCF="T" alternateAlleleVCF="{alt}"/>
      </Location></SimpleAllele>
      <Classifications><GermlineClassification DateLastEvaluated="2025-01-02"><Description>{classification}</Description><ReviewStatus>criteria provided, conflicting classifications</ReviewStatus></GermlineClassification></Classifications>
      <RCVList><RCVAccession Accession="RCV000000001" Version="3"><ClassifiedConditionList><ClassifiedCondition>Test condition</ClassifiedCondition></ClassifiedConditionList><RCVClassifications><GermlineClassification><Description DateLastEvaluated="2025-01-01">Uncertain significance</Description><ReviewStatus>criteria provided, single submitter</ReviewStatus></GermlineClassification></RCVClassifications></RCVAccession></RCVList>
    </ClassifiedRecord></VariationArchive>''')


class ProjectionTest(unittest.TestCase):
    def setUp(self):
        self.target = SimpleNamespace(gene='INS')
        self.mock = {'location': {'start': 1, 'end': 6, 'strand': 1}, 'sequence': 'ATGTAA',
                     'protein': {'segments': [{'start': 1, 'end': 6}], 'translation': 'M'}}
        self.impact = {'chromosome': 'chr11', 'complemented': True,
                       'runs': [{'local': 1, 'genomic': 102, 'step': -1, 'length': 6}]}
        self.root = ET.Element('ClinVarResult-Set')
        self.root.append(record())

    def project(self):
        return project_records(self.root, self.target, self.mock, self.impact)

    def test_reverse_genomic_coordinates_and_complement_are_independent(self):
        variants, excluded = self.project()
        self.assertFalse(excluded)
        v = variants[0]
        self.assertEqual((v['position'], v['ref'], v['alt']), (1, 'A', 'G'))
        self.assertEqual((v['coding_change'], v['protein_change']), ('c.1A>G', 'p.M1V'))
        self.assertEqual(v['consequence'], 'start codon')
        self.assertEqual(v['classification'], 'Conflicting classifications of pathogenicity')
        self.assertEqual(v['conditions'][0]['classification'], 'Uncertain significance')

    def test_negative_record_sequence_is_indexed_from_end_without_second_complement(self):
        self.mock['location']['strand'] = -1
        self.impact['runs'] = [{'local': 1, 'genomic': 97, 'step': 1, 'length': 6}]
        variants, _ = self.project()
        self.assertEqual((variants[0]['position'], variants[0]['ref'], variants[0]['protein_change']), (6, 'A', 'p.M1V'))

    def test_two_alternatives_at_one_position_remain_two_records(self):
        self.root.append(record('2', 'G', 'Benign'))
        variants, _ = self.project()
        self.assertEqual({v['alt'] for v in variants}, {'C', 'G'})
        self.assertEqual(len(variants), 2)

    def test_gene_span_cannot_supply_allele_location(self):
        self.root[0].find('ClassifiedRecord/SimpleAllele/Location').clear()
        variants, excluded = self.project()
        self.assertFalse(variants)
        self.assertEqual(excluded, {'no_unique_GRCh38_location': 1})

    def test_reference_mismatch_is_excluded_not_reassigned(self):
        self.mock['sequence'] = 'GTGTAA'
        self.mock['protein']['translation'] = 'V'
        self.assertEqual(self.project()[1], {'reference_mismatch': 1})

    def test_indels_and_missing_germline_classification_are_scoped_out(self):
        self.root[0].find('ClassifiedRecord/SimpleAllele/VariantType').text = 'Deletion'
        self.assertEqual(self.project()[1], {'not_a_single_base_substitution': 1})
        self.root[0] = record()
        self.root[0].find('ClassifiedRecord/Classifications').clear()
        self.assertEqual(self.project()[1], {'no_germline_classification': 1})

    def test_a_shortened_intron_middle_is_not_outside_the_gene(self):
        # Two runs with the chromosome jumping between them, as a gene drawn
        # with shortened introns is mapped: 102-100, then 90-88.
        self.mock['sequence'] = 'ATGAAATAA'
        self.mock['location']['end'] = 9
        self.mock['protein'] = {'segments': [{'start': 1, 'end': 9}], 'translation': 'MK'}
        self.impact['runs'] = [{'local': 1, 'genomic': 102, 'step': -1, 'length': 3},
                               {'local': 4, 'genomic': 90, 'step': -1, 'length': 6}]
        self.root[0] = record('1', at=95)
        self.root.append(record('2', at=150))
        self.root.append(record('3', at=80))
        variants, excluded = self.project()
        self.assertFalse(variants)
        self.assertEqual(excluded, {'intron_not_drawn': 1, 'outside_drawn_gene': 2})

    def test_one_run_leaves_no_intron_undrawn(self):
        self.root[0] = record('1', at=200)
        self.assertEqual(self.project()[1], {'outside_drawn_gene': 1})

    def test_duplicate_variation_id_aborts(self):
        self.root.append(copy.deepcopy(self.root[0]))
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            self.project()

    def test_replay_requires_every_searched_id_and_intact_source_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            raw = ET.tostring(self.root)
            (cache / 'batch-0.xml').write_bytes(raw)
            metadata = {'query': 'INS[gene]', 'record_count': 2, 'ids': ['1', '2'],
                        'batches': [{'file': 'batch-0.xml', 'sha256': hashlib.sha256(raw).hexdigest()}]}
            (cache / 'source.json').write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError, 'all search IDs'):
                fetch('INS', cache, replay=True)
            metadata.update(record_count=1, ids=['1'])
            (cache / 'source.json').write_text(json.dumps(metadata))
            self.assertEqual(len(fetch('INS', cache, replay=True)[0]), 1)
            (cache / 'batch-0.xml').write_bytes(raw + b' ')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                fetch('INS', cache, replay=True)


class Response:
    """What `Entrez.efetch` hands back: a context manager with a body."""

    def __init__(self, body):
        self.body = body

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def read(self):
        if isinstance(self.body, Exception):
            raise self.body
        return self.body


def batch_of(*ids):
    return ('<ClinVarResult-Set>' + ''.join(
        f'<VariationArchive VariationID="{i}"/>' for i in ids) + '</ClinVarResult-Set>').encode()


class BatchTest(unittest.TestCase):
    def efetch(self, *bodies):
        calls = []

        def efetch(**kwargs):
            calls.append(kwargs['id'])
            return Response(bodies[len(calls) - 1])
        return efetch, calls

    def test_a_dropped_read_and_a_short_batch_are_fetched_again(self):
        efetch, calls = self.efetch(http.client.IncompleteRead(b'<Clin'), batch_of('1'), batch_of('2', '1'))
        waits = []
        raw = fetch_batch(efetch, ['1', '2'], pause=waits.append)
        self.assertEqual(raw, batch_of('2', '1'))
        self.assertEqual(calls, ['1,2', '1,2', '1,2'])
        self.assertEqual(waits, [10, 20])

    def test_a_batch_that_never_arrives_whole_stops_the_bake(self):
        efetch, calls = self.efetch(batch_of('1'), batch_of('1'))
        with self.assertRaisesRegex(ValueError, '1 of the 2 records'):
            fetch_batch(efetch, ['1', '2'], attempts=2, pause=lambda _: None)
        self.assertEqual(len(calls), 2)

    def test_an_error_response_is_not_a_batch(self):
        efetch, _ = self.efetch(b'<ClinVarResult-Set><ERROR>busy</ERROR></ClinVarResult-Set>', batch_of('1'))
        self.assertEqual(fetch_batch(efetch, ['1'], pause=lambda _: None), batch_of('1'))


class ChosenTest(unittest.TestCase):
    def test_all_means_every_row_with_a_snapshot(self):
        every = chosen(SimpleNamespace(all=True, target=None))
        self.assertTrue(every)
        self.assertTrue(all(t.clinvar_available for t in every))
        self.assertEqual([t.slug for t in chosen(SimpleNamespace(all=False, target='insulin'))], ['insulin'])


def traited(id='1', medgen='C3150617', rcv_names=(('Maturity-onset diabetes of the young type 10', 'C3150617'),)):
    """A mapped record whose RCV names its conditions and whose ConditionList
    carries the identifiers ClinVar gives them."""
    conditions = ''.join(
        f'<ClassifiedCondition DB="MedGen" ID="{cui}">{name}</ClassifiedCondition>' if cui
        else f'<ClassifiedCondition>{name}</ClassifiedCondition>'
        for name, cui in rcv_names)
    return ET.fromstring(f'''<VariationArchive VariationID="{id}" Accession="VCV00000000{id}" Version="1" VariationName="synthetic">
    <ClassifiedRecord>
      <RCVList><RCVAccession Accession="RCV00000000{id}" Version="1">
        <ClassifiedConditionList>{conditions}</ClassifiedConditionList>
        <RCVClassifications><GermlineClassification><Description>Likely pathogenic</Description></GermlineClassification></RCVClassifications>
      </RCVAccession></RCVList>
      <Classifications><GermlineClassification><ConditionList>
        <TraitSet Type="Disease"><Trait Type="Disease">
          <Name><ElementValue Type="Preferred">Maturity-onset diabetes of the young type 10</ElementValue><XRef ID="MONDO:0013240" DB="MONDO"/></Name>
          <Symbol><ElementValue Type="Preferred">MODYX</ElementValue></Symbol>
          <Symbol><ElementValue Type="Alternate">MODY10</ElementValue><XRef Type="MIM" ID="613370" DB="OMIM"/></Symbol>
          <XRef ID="{medgen}" DB="MedGen"/><XRef ID="MONDO:0013240" DB="MONDO"/><XRef Type="MIM" ID="613370" DB="OMIM"/>
        </Trait></TraitSet>
        <TraitSet Type="Disease"><Trait Type="Disease">
          <Name><ElementValue Type="Preferred">Permanent neonatal diabetes mellitus</ElementValue></Name>
          <Symbol><ElementValue Type="Alternate">PDMI</ElementValue><XRef Type="MIM" ID="606176" DB="OMIM"/></Symbol>
          <Symbol><ElementValue Type="Preferred">PNDM</ElementValue></Symbol>
          <XRef ID="C1833104" DB="MedGen"/><XRef Type="Phenotypic series" ID="PS606176" DB="OMIM"/>
        </Trait></TraitSet>
      </ConditionList></GermlineClassification></Classifications>
    </ClassifiedRecord></VariationArchive>''')


class TraitTest(unittest.TestCase):
    def test_a_condition_takes_its_own_entry_symbol_and_identifiers(self):
        root = ET.Element('ClinVarResult-Set')
        root.append(traited(rcv_names=(
            ('Maturity-onset diabetes of the young type 10', 'C3150617'),
            ('Permanent neonatal diabetes mellitus', 'C1833104'),
            ('INS-related disorder', None),
            ('not provided', 'C3661900'),
        )))
        self.assertEqual(condition_traits(root, [{'variation_id': '1'}]), {
            # The symbol OMIM cross-references for the entry, not the preferred one.
            'Maturity-onset diabetes of the young type 10': {
                'medgen': 'C3150617', 'symbol': 'MODY10', 'omim': '613370', 'mondo': 'MONDO:0013240'},
            # A phenotypic series names no symbol of its own: the preferred one stands.
            'Permanent neonatal diabetes mellitus': {
                'medgen': 'C1833104', 'symbol': 'PNDM', 'omim': 'PS606176'},
        })

    def test_unmapped_records_contribute_nothing(self):
        root = ET.Element('ClinVarResult-Set')
        root.append(traited())
        self.assertEqual(condition_traits(root, [{'variation_id': '2'}]), {})

    def test_symbol_order_in_a_record_does_not_decide(self):
        def t2d(id, order):
            symbols = {
                'T2D': '<Symbol><ElementValue Type="Alternate">T2D</ElementValue><XRef Type="MIM" ID="125853" DB="OMIM"/></Symbol>',
                'NIDDM': '<Symbol><ElementValue Type="Alternate">NIDDM</ElementValue><XRef ID="HP:0005978" DB="Human Phenotype Ontology"/><XRef Type="MIM" ID="125853" DB="OMIM"/></Symbol>',
            }
            return ET.fromstring(f'''<VariationArchive VariationID="{id}"><ClassifiedRecord>
              <RCVList><RCVAccession><ClassifiedConditionList><ClassifiedCondition DB="MedGen" ID="C0011860">Type 2 diabetes mellitus</ClassifiedCondition></ClassifiedConditionList>
                <RCVClassifications><GermlineClassification/></RCVClassifications></RCVAccession></RCVList>
              <Classifications><GermlineClassification><ConditionList><TraitSet><Trait>
                {''.join(symbols[s] for s in order)}
                <XRef ID="C0011860" DB="MedGen"/><XRef Type="MIM" ID="125853" DB="OMIM"/>
              </Trait></TraitSet></ConditionList></GermlineClassification></Classifications>
            </ClassifiedRecord></VariationArchive>''')
        root = ET.Element('ClinVarResult-Set')
        root.append(t2d('1', ('NIDDM', 'T2D')))
        root.append(t2d('2', ('T2D', 'NIDDM')))
        self.assertEqual(
            condition_traits(root, [{'variation_id': '1'}, {'variation_id': '2'}])['Type 2 diabetes mellitus']['symbol'],
            'T2D')

    def test_one_name_with_two_concepts_aborts(self):
        root = ET.Element('ClinVarResult-Set')
        root.append(traited('1'))
        root.append(traited('2', medgen='C0000001', rcv_names=(
            ('Maturity-onset diabetes of the young type 10', 'C0000001'),)))
        with self.assertRaisesRegex(ValueError, 'conflicting identifiers'):
            condition_traits(root, [{'variation_id': '1'}, {'variation_id': '2'}])


if __name__ == '__main__':
    unittest.main()
