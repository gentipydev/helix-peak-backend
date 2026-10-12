"""The bake of an experimental entry chosen by rule: model, scene and fold together.

Offline: RCSB, PDBe, UniProt and InterPro answer from `fixtures/pdb/`. What
needs no tool runs anywhere numpy is. The whole bake also needs scipy and
trimesh, PyMOL, the Dart SDK and the app's checkout for its scene importer,
and the app's check needs Flutter; each skips without what it needs, as
`test_alphafold.py`'s do.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import shutil
import stat
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

np = pytest.importorskip("numpy")

from pipeline.paths import CLIENT  # noqa: E402
from pipeline.structure import (  # noqa: E402
    alphafold, experimental, fold_job, normalize, pdb_policy, pdb_search)
from pipeline.structure.normalize import Unsuitable  # noqa: E402
from pipeline.structure.test_fold_job import NOW, services  # noqa: E402
from pipeline.targets import BY_SLUG, Region, Source, Target  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "pdb"
RLX = (FIXTURES / "6rlx_cut.cif").read_bytes()
MV1 = (FIXTURES / "2mv1_cut.cif").read_bytes()

RELAXIN = BY_SLUG["relaxin"]
UNIPROT = json.loads((FIXTURES / "uniprot_P04090.json").read_bytes())
# The record as the bakes read it: the gene, and the protein it makes.
RECORD = json.dumps({"gene": "RLN2", "protein": {"translation": UNIPROT["sequence"]["value"]},
                     "peptides": []}).encode()
ROW = {"slug": "relaxin", "display": "Relaxin", "gene": "RLN2", "uniprot": "P04090",
       "accession": "NC_000009.12", "summary": "A hormone of two chains.",
       "facts": {"residues": 185, "exons": 2, "chains": 2, "bridges": 3},
       "catalog_order": None, "chain": None, "chains": [], "structure": None, "tracks": {}}


def relaxin():
    """Relaxin's inputs, its best candidate (6RLX, both chains), and that
    candidate held to its file."""
    inputs = fold_job.inputs_of(RELAXIN, UNIPROT["sequence"]["value"], UNIPROT)
    whole = json.loads((FIXTURES / "graphql_relaxin.json").read_bytes())["data"]
    entities = [pdb_search._entity(e) for e in whole["polymer_entities"]]
    candidate = pdb_policy.rank(inputs, entities).candidates[0]
    entry = next(e.entry for e in entities if e.entry.entry == candidate.entry)
    return inputs, candidate, normalize.examine(inputs, candidate, RLX), entry


# ------------------------------------------------------------ what needs no tool


def test_an_authors_numbering_is_an_offset_only_where_one_number_does_for_every_residue():
    def chain(*numbered):
        return SimpleNamespace(residues=[SimpleNamespace(number=n, auth=a) for n, a in numbered])

    assert experimental.author_offset(chain((25, "1"), (26, "2"), (27, "3"))) == 24
    assert experimental.author_offset(chain((25, "-3"), (26, "-2"))) == 28
    # An insertion code, and a jump in the author's count, are no offset.
    assert experimental.author_offset(chain((25, "1"), (26, "1A"), (27, "2"))) is None
    assert experimental.author_offset(chain((25, "1"), (26, "2"), (27, "10"))) is None


def test_relaxins_chains_are_numbered_by_the_entry_from_one_and_the_track_says_by_how_much():
    _, _, examined, _ = relaxin()
    b_chain, a_chain = examined.chains
    assert (experimental.author_offset(b_chain), experimental.author_offset(a_chain)) == (28, 165)


def test_a_target_crosses_to_the_bake_as_json_and_comes_back_itself():
    said = json.loads(json.dumps(experimental.target_json(RELAXIN)))
    assert said["structure"] is None
    assert experimental.target_of_json(said) == dataclasses.replace(RELAXIN, structure=None)


def _built_row() -> dict:
    """A protein row as the resolver writes one (`resolve.resolve`)."""
    return {
        "slug": "ttr", "gene": "TTR", "uniprot": "P02766", "accession": "NG_009490.1",
        "slice_start": None, "slice_end": None, "slice_strand": None,
        "transcript_id": "NM_000371.4", "protein_id": "NP_000362.1", "taxon_id": 9606,
        "display": "Transthyretin", "summary": "Carries thyroxine.", "chain_name": None,
        "mature_peptides": False, "residues": 147, "exons": 4, "chains": 1, "bridges": 0,
        "regions": [], "disulfides": [], "structure": None, "resolver_version": 1,
        "catalog_order": None,
        "provenance": {
            "cleaved": True, "chain_label": "Transthyretin", "uniprot_variants": [[26, "S", "G"]],
            "regions_from": [
                {"label": "Signal peptide", "short": "S", "start": 1, "end": 20, "kept": False},
                {"label": "Transthyretin", "short": "T", "start": 21, "end": 147, "kept": True}],
        },
    }


def test_a_row_is_handed_to_the_apps_check_as_the_service_serves_it():
    # The service reads its settings as it is imported (`tests/conftest.py`).
    os.environ.setdefault("NCBI_EMAIL", "tests@example.com")
    from app import catalog

    row = _built_row()
    row["structure"] = {"chrome": {"pdb": "AF-P02766-F1"}, "chains": [{"node": "x", "tint": "y"}]}
    columns = [c.strip().split(".")[-1] for c in catalog._SUMMARY_COLUMNS.split(",")]
    served = catalog._summary(tuple(row[column] for column in columns))
    assert experimental.served_row(row) == served
    # And without a fold yet, as a protein just resolved is.
    assert experimental.served_row(_built_row())["structure"] is None
    assert experimental.served_row(_built_row())["chains"] == []


def test_a_built_proteins_target_is_rebuilt_from_what_the_service_serves():
    from pipeline.resolver.resolve import target_of

    row = _built_row()
    served = {**experimental.served_row(row),
              **{key: row[key] for key in ("mature_peptides", "transcript_id", "protein_id",
                                           "regions", "disulfides", "provenance")}}
    assert experimental.target_of_served(served) == target_of(row)
    # One of the twenty has no such table on its row, and is not baked this way.
    with pytest.raises(ValueError, match="not a protein built on demand"):
        experimental.target_of_served({**served, "provenance": {"cleaved": True}})


def test_the_bake_is_handed_the_proteins_own_row_with_the_chains_as_its_structure():
    inputs, candidate, examined, _ = relaxin()
    drawn = experimental._bake_target(RELAXIN, inputs, candidate, examined, 8, "cartoon")
    assert drawn.structure.pdb == "6RLX" and drawn.structure.bonds
    assert [(c.node, c.pdb_chain) for c in drawn.structure.chains] == [
        ("chainA", "A"), ("chainB", "B")]
    assert drawn.structure.residues is None and drawn.structure.sampling == 8
    assert drawn.regions == RELAXIN.regions and drawn.slug == "relaxin"


def test_a_cut_protein_with_no_kept_region_is_baked_as_what_its_removed_ones_leave():
    # Apelin's row names what is removed and nothing as kept; the products
    # are what is left, and the bakers read a cut precursor by its kept regions.
    target = Target(
        slug="apln", gene="APLN", uniprot="Q9ULZ1", display="Apelin", source=Source("NG_1.1"),
        structure=None, cleaved=True, chain_label="Apelin", aa=30,
        regions=(Region("Signal peptide", "S", 1, 10, kept=False),))
    inputs, candidate, examined, _ = relaxin()
    inputs = dataclasses.replace(
        inputs, products=(pdb_policy.Product(0, "Apelin", 11, 30),), translation="M" * 30)
    drawn = experimental._bake_target(target, inputs, candidate, examined, 8, "cartoon")
    assert [(r.label, r.start, r.end, r.kept) for r in drawn.regions] == [
        ("Signal peptide", 1, 10, False), ("Apelin", 11, 30, True)]


# ------------------------------------------------------------ the app's check


def _flutter(tmp_path: Path, verdict, status: int = 0, sleep: int = 0) -> str:
    """A stand-in for `flutter test`: it writes `verdict` (None for nothing)
    beside the bundle it is handed, and exits with `status`."""
    script = tmp_path / "flutter"
    script.write_text(
        "#!" + sys.executable + "\n"
        "import json, os, sys, time\n"
        f"time.sleep({sleep})\n"
        "assert sys.argv[1:] == ['test', '--no-pub', 'test/fold_bundle_check.dart'], sys.argv\n"
        f"verdict = {verdict!r}\n"
        "if verdict is not None:\n"
        "    open(os.path.join(os.environ['FOLD_BUNDLE'], 'check.json'), 'w').write(\n"
        "        verdict if isinstance(verdict, str) else json.dumps(verdict))\n"
        f"sys.exit({status})\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return str(script)


def _folder(tmp_path: Path) -> Path:
    folder = tmp_path / "bundle"
    folder.mkdir(exist_ok=True)
    return folder


def test_a_fold_the_app_holds_is_returned_with_what_the_check_measured(tmp_path):
    measured = {"verdict": "ok", "steps": 6, "rings": 406, "bind_ms": 48}
    assert experimental.app_check(_folder(tmp_path), tmp_path,
                                  _flutter(tmp_path, measured)) == measured


def test_a_rule_of_the_animation_the_entry_breaks_is_a_refusal_in_the_checks_own_words(tmp_path):
    said = {"verdict": "unsuitable", "steps": 2,
            "failed": {"step": "every ring of its model lies on its own stretch of chain",
                       "detail": "A ring of chainA lies on\n  no stretch of its chain"}}
    with pytest.raises(Unsuitable) as refused:
        experimental.app_check(_folder(tmp_path), tmp_path, _flutter(tmp_path, said, status=1))
    assert refused.value.reason == "binding"
    assert refused.value.detail == ("every ring of its model lies on its own stretch of chain: "
                                    "A ring of chainA lies on no stretch of its chain")


@pytest.mark.parametrize("verdict, status, said", [
    # The bundle is not one the app can read: a fault of the bake's, or the check's.
    ({"verdict": "fault", "failed": {"step": "(reading the bundle)", "detail": "bad scene"}},
     1, "broke.*bad scene"),
    # It ended without saying how.
    (None, 1, "gave no verdict"),
    ("not json", 1, "gave no verdict"),
    # What it says and how it ended disagree.
    ({"verdict": "ok", "steps": 6}, 1, "broke"),
    ({"verdict": "unsuitable", "failed": {"step": "s", "detail": "d"}}, 0, "broke"),
])
def test_a_check_that_breaks_is_no_verdict_on_the_entry(tmp_path, verdict, status, said):
    with pytest.raises(RuntimeError, match=said):
        experimental.app_check(_folder(tmp_path), tmp_path, _flutter(tmp_path, verdict, status))


def test_a_verdict_left_by_an_earlier_run_is_never_read(tmp_path):
    folder = _folder(tmp_path)
    (folder / "check.json").write_text(json.dumps({"verdict": "ok", "steps": 6}))
    with pytest.raises(RuntimeError, match="gave no verdict"):
        experimental.app_check(folder, tmp_path, _flutter(tmp_path, None))


def test_a_check_that_does_not_end_or_cannot_start_is_no_verdict_either(tmp_path):
    with pytest.raises(RuntimeError, match="had not finished after 1 s"):
        experimental.app_check(_folder(tmp_path), tmp_path,
                               _flutter(tmp_path, {"verdict": "ok"}, sleep=5), timeout=1)
    with pytest.raises(RuntimeError, match="could not be started"):
        experimental.app_check(_folder(tmp_path), tmp_path, str(tmp_path / "no-flutter-here"))


def test_the_bundle_is_the_row_with_this_folds_chains_and_words(tmp_path):
    built = experimental.Built(
        glb=b"glb", scene=b"scene", folding=b"{}", chrome={"pdb": "6RLX"},
        chains=[{"node": "chainA", "tint": "mature1"}], provenance={}, folding_provenance={},
        working_copy="", audit={})
    stale = {**ROW, "chains": [{"node": "plddtVeryHigh", "tint": "plddtVeryHigh"}],
             "structure": {"pdb": "AF-P04090-F1"}}
    folder = experimental.bundle(built, stale, tmp_path / "b")
    row = json.loads((folder / "row.json").read_text())
    assert row["structure"] == {"pdb": "6RLX"} and row["chains"] == built.chains
    assert row["gene"] == "RLN2"
    assert (folder / "model.glb").read_bytes() == b"glb"
    assert (folder / "model.fsceneb").read_bytes() == b"scene"
    assert (folder / "folding.json").read_bytes() == b"{}"


# ------------------------------------------------------------ the whole bake


def _bake_tools() -> str:
    """PyMOL, Dart and the app's checkout, or why the whole bake cannot run here."""
    pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    from pipeline.structure import bake

    if not Path(bake.PYMOL).exists():
        pytest.skip("PyMOL is not installed")
    sdk = Path.home() / "flutter" / "bin" / "cache" / "dart-sdk" / "bin" / "dart"
    dart = str(sdk) if sdk.exists() else shutil.which("dart")
    if not dart:
        pytest.skip("no Dart SDK to run the scene importer with")
    unready = alphafold.importer_unready(CLIENT)
    if unready:
        pytest.skip(unready)
    return dart


def _flutter_tool() -> str:
    flutter = Path.home() / "flutter" / "bin" / "flutter"
    found = str(flutter) if flutter.exists() else shutil.which("flutter")
    if not found:
        pytest.skip("no Flutter to run the app's check with")
    if not (CLIENT / experimental.BUNDLE_CHECK).exists():
        pytest.skip(f"the app's checkout has no {experimental.BUNDLE_CHECK}")
    return found


@pytest.fixture(scope="module")
def baked(tmp_path_factory):
    """Relaxin, baked once from 6RLX's fixture."""
    dart = _bake_tools()
    inputs, candidate, examined, entry = relaxin()
    built = experimental.build(RELAXIN, inputs, candidate, examined, RECORD, entry, None,
                               tmp_path_factory.mktemp("relaxin"), CLIENT, dart)
    return SimpleNamespace(built=built, inputs=inputs, candidate=candidate, examined=examined,
                           entry=entry, dart=dart)


def test_relaxin_is_baked_as_the_twenty_are_model_and_fold_from_one_export(baked):
    from pipeline.structure.glb import read_glb

    built = baked.built
    meshes = read_glb(built.glb)
    assert list(meshes) == ["chainA", "chainB", "bonds"]
    assert built.chains == [{"node": "chainA", "tint": "mature1"},
                            {"node": "chainB", "tint": "mature2"},
                            {"node": "bonds", "tint": "cysteine"}]
    everything = np.vstack([m.positions for m in meshes.values()])
    assert abs(np.ptp(everything, axis=0).max() - 1.0) < 1e-6      # cut to the twenty's measure
    assert 0 < len(built.scene) <= alphafold.BUDGET_BYTES

    track = json.loads(built.folding)
    # One snapshot and one frame: the track names the very model it was made with.
    assert track["frame"]["glb_sha256"] == hashlib.sha256(built.glb).hexdigest()
    assert track["frame"]["matched_within"] == 0.0
    assert (track["slug"], track["gene"], track["uniprot"], track["pdb"]) == (
        "relaxin", "RLN2", "P04090", "6RLX")
    assert track["schema_version"] == 2 and len(track["bridges"]) == 3
    b_chain, a_chain = track["chains"]
    # Each chain as the entry names it, and its author's numbering beside the precursor's.
    assert (b_chain["node"], b_chain["pdb_chain"], b_chain["offset"],
            b_chain["first"], b_chain["last"]) == ("chainA", "B", 28, 25, 53)
    assert (a_chain["node"], a_chain["pdb_chain"], a_chain["offset"],
            a_chain["first"], a_chain["last"]) == ("chainB", "A", 165, 162, 185)
    # The whole product is listed. The construct began at Ser26, so Asp25 was
    # never in the experiment; and its last residue, Ser53, has no place in it.
    states = [r["state"] for r in b_chain["residues"]]
    assert states[0] == "absent" and states[-1] == "disordered"
    assert states.count("ordered") == 27 and len(states) == 29
    assert all(r["state"] == "ordered" for r in a_chain["residues"])
    assert track["secondary_structure"].startswith("DSSP in the entry's file from PDBe")
    assert track["built_by"] == "pipeline/structure/experimental.py"


def test_the_fold_pages_words_and_the_provenance_are_the_entrys_own_numbers(baked):
    built = baked.built
    assert built.chrome["pdb"] == "6RLX" and built.chrome["label"] == "the mature chains"
    assert built.chrome["count"] == 53 and built.chrome["modelled"] is None
    assert built.chrome["sentence"] == "X-ray, 1.50 Å. Two chains, joined by two disulfide bridges."
    said = built.provenance
    assert (said["pdb"], said["source"], said["licence"], said["policy"]) == (
        "6RLX", "PDB", "CC0 1.0", pdb_policy.POLICY_VERSION)
    assert said["nodes"] == ["chainA", "chainB", "bonds"]
    assert said["method"] == ["X-RAY DIFFRACTION"] and said["resolution"] == 1.5
    assert said["file"]["sha256"] == hashlib.sha256(RLX).hexdigest()
    assert said["residues"] == {"ordered": 51, "disordered": 1, "absent": 1}
    assert said["bridges"] == [[35, 172], [47, 185], [171, 176]]
    assert said["pymol"] == experimental.PYMOL and said["representation"] == "cartoon"
    assert said["importer"] == {"package": "flutter_scene", "version": alphafold.FLUTTER_SCENE}
    assert [c["auth"] for c in said["chains"]] == ["B", "A"]
    assert all(frame["ca_max"] <= alphafold.FRAME_NEAR for frame in said["frame"].values())
    assert built.folding_provenance["frame"]["glb_sha256"] == hashlib.sha256(built.glb).hexdigest()
    # The audit names what was read, and keeps the whole map.
    audit = built.audit
    assert audit["file_sha256"] == said["file"]["sha256"]
    assert audit["record_sha256"] == hashlib.sha256(RECORD).hexdigest()
    # 6RLX numbers its B chain from -2, and has no place for its last residue.
    assert len(audit["map"]) == 52 and audit["map"][0] == {
        "node": "chainA", "n": 26, "label": 1, "auth": "-2", "residue": "SER", "letter": "S",
        "placed": True}
    assert audit["map"][27] == {
        "node": "chainA", "n": 53, "label": 28, "auth": "25", "residue": "SER", "letter": "S",
        "placed": False}
    assert audit["working_copy_sha256"] == hashlib.sha256(built.working_copy.encode()).hexdigest()


def test_the_same_entry_bakes_the_same_bytes(baked, tmp_path):
    again = experimental.build(RELAXIN, baked.inputs, baked.candidate, baked.examined, RECORD,
                               baked.entry, None, tmp_path, CLIENT, baked.dart)
    assert (again.glb, again.scene, again.folding) == (
        baked.built.glb, baked.built.scene, baked.built.folding)


def test_a_model_no_sampling_fits_in_the_budget_is_refused_for_its_size(baked, tmp_path):
    with pytest.raises(Unsuitable) as refused:
        experimental.build(RELAXIN, baked.inputs, baked.candidate, baked.examined, RECORD,
                           baked.entry, None, tmp_path, CLIENT, baked.dart, budget=1000)
    assert refused.value.reason == "size"


def test_an_entry_that_assigns_no_shape_is_given_pymols_and_the_track_shares_it(baked, tmp_path):
    bare = dataclasses.replace(baked.examined, secondary={"chainA": (), "chainB": ()},
                               secondary_method=None)
    built = experimental.build(RELAXIN, baked.inputs, baked.candidate, bare, RECORD, baked.entry,
                               None, tmp_path, CLIENT, baked.dart)
    assert built.provenance["secondary_structure"] == f"PyMOL {experimental.PYMOL} dss"
    track = json.loads(built.folding)
    assert track["secondary_structure"] == f"PyMOL {experimental.PYMOL} dss, 6RLX"
    shapes = {r["ss"] for c in track["chains"] for r in c["residues"] if r["state"] == "ordered"}
    assert "helix" in shapes
    # What PyMOL assigned is what the copy both bakes read says.
    assert "HELIX" in built.working_copy and built.provenance["representation"] == "cartoon"


def test_a_ribbon_measured_where_the_check_would_not_have_it_is_not_kept(baked):
    track = json.loads(baked.built.folding)
    measured = [r for c in track["chains"] for r in c["residues"] if "ribbon" in r]
    assert measured and experimental._unmeasured(track, baked.built.glb) == 0
    # A centre moved off its ribbon, as a neighbouring strand's would be.
    moved = measured[0]
    moved["ribbon"] = [moved["ribbon"][0] + 0.3] + moved["ribbon"][1:]
    assert experimental._unmeasured(track, baked.built.glb) == 1
    assert "ribbon" not in moved and "ribbon" in measured[1]


def test_the_baked_fold_holds_to_the_apps_own_rules(baked, tmp_path):
    flutter = _flutter_tool()
    folder = experimental.bundle(baked.built, ROW, tmp_path / "bundle")
    measured = experimental.app_check(folder, CLIENT, flutter)
    assert measured["verdict"] == "ok" and measured["steps"] == 6
    assert measured["residues"] == 52 and measured["bridges"] == 3
    assert measured["largest_nudge"] < 0.5 and measured["largest_step"] < 1.1


def test_a_fold_the_binder_cannot_lay_on_its_chain_is_refused_by_the_app(baked, tmp_path):
    flutter = _flutter_tool()
    track = json.loads(baked.built.folding)
    for residue in track["chains"][0]["residues"]:
        if "ca" in residue:
            residue["ca"] = [residue["ca"][0] + 0.6] + residue["ca"][1:]
    moved = dataclasses.replace(baked.built, folding=json.dumps(track).encode())
    with pytest.raises(Unsuitable) as refused:
        experimental.app_check(experimental.bundle(moved, ROW, tmp_path / "bundle"), CLIENT, flutter)
    assert refused.value.reason == "binding"
    assert "lies on no stretch of its chain" in refused.value.detail


def test_a_protein_is_searched_chosen_baked_and_checked_in_one_go(tmp_path):
    dart, flutter = _bake_tools(), _flutter_tool()
    post, get, asked = services(files={"6RLX": RLX, "2MV1": MV1})
    made = experimental.make(
        RELAXIN, RECORD, ROW, tmp_path / "work", CLIENT, dart, flutter,
        fetch_entry=lambda accession: SimpleNamespace(body=UNIPROT, release="2026_03"),
        pfam=lambda accession: (), get=get, post=post, now=lambda: NOW)
    assert made.choice.outcome == "chosen" and made.uniprot_release == "2026_03"
    assert fold_job.candidate_id(made.choice.winner.candidate) == "6RLX:B@26,A@162"
    assert asked.files == ["6RLX"]
    built = made.built
    assert built.measured["scene_bytes"] == len(built.scene) and built.measured["rings"] > 0

    said = experimental.write(made, tmp_path / "out")
    out = tmp_path / "out"
    assert json.loads((out / "outcome.json").read_text()) == said
    assert said["outcome"] == "chosen" and said["candidate"] == "6RLX:B@26,A@162"
    assert said["sha256"] == {
        "glb": hashlib.sha256((out / "model.glb").read_bytes()).hexdigest(),
        "scene": hashlib.sha256((out / "model.fsceneb").read_bytes()).hexdigest(),
        "folding": hashlib.sha256((out / "folding.json").read_bytes()).hexdigest()}
    assert said["folding_provenance"]["frame"]["glb_sha256"] == said["sha256"]["glb"]
    audit = json.loads((out / "audit.json").read_text())
    assert audit["selection"]["outcome"] == "chosen" and audit["selection"]["policy"] == 1
    assert audit["selection"]["search"]["entities"] == 4
    assert audit["measured"] == built.measured


def test_a_search_that_finds_nothing_leaves_no_fold_and_says_so(tmp_path):
    post, get, _ = services(ids=[])
    made = experimental.make(
        RELAXIN, RECORD, ROW, tmp_path / "work", CLIENT, "dart", "flutter",
        fetch_entry=lambda accession: SimpleNamespace(body=UNIPROT, release="2026_03"),
        pfam=lambda accession: (), get=get, post=post, now=lambda: NOW, recovery=False)
    assert made.choice.outcome == "none" and made.built is None
    said = experimental.write(made, tmp_path / "out")
    assert said["outcome"] == "none" and "chrome" not in said
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == ["outcome.json"]


def test_an_outage_while_a_protein_is_made_is_raised(tmp_path):
    def down(accession):
        raise OSError("UniProt is not answering")

    with pytest.raises(OSError):
        experimental.make(RELAXIN, RECORD, ROW, tmp_path, CLIENT, "dart", "flutter",
                          fetch_entry=down, pfam=lambda accession: ())
    assert not os.listdir(tmp_path)
