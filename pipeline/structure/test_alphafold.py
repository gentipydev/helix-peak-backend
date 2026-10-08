"""An AlphaFold DB model as a protein's fold: what is taken, drawn and refused.

The network is never reached: the API's answers and one real model are under
`fixtures/alphafold/`. The rules need numpy alone and run with the service's
suite. The mesh tests need the structure bake's environment (trimesh, scipy)
and skip without it; the whole bake also needs PyMOL and the app's checkout
for its scene importer, and skips without those.
"""

from __future__ import annotations

import json
import shutil
import sys
import urllib.error
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import CLIENT  # noqa: E402
from pipeline.structure import alphafold  # noqa: E402
from pipeline.structure.alphafold import (  # noqa: E402
    BANDS, Confidence, Entry, Job, Model, Refused)
from pipeline.structure.pdb import ca_atoms, ssbonds  # noqa: E402
from pipeline.targets import Chain, Structure  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "alphafold"
SLN = "MGINTRELFLNFTIVLITVILMWLLVRSYQY"
SLN_MODEL = FIXTURES / "AF-O00631-F1-model_v6.pdb"
SLN_ENTRY = Entry("O00631", "AF-O00631-F1", 6,
                  "https://alphafold.ebi.ac.uk/files/AF-O00631-F1-model_v6.pdb")


def _answer(accession: str) -> list:
    return json.loads((FIXTURES / f"{accession}.json").read_text())


def _serving(pages: dict):
    """A `get` that answers from `pages` and fails on any other address."""
    def get(url: str) -> bytes:
        if url not in pages:
            raise AssertionError(f"the network was asked for {url}")
        page = pages[url]
        if isinstance(page, Exception):
            raise page
        return page
    return get


def _not_found(url: str) -> urllib.error.HTTPError:
    return urllib.error.HTTPError(url, 404, "Not Found", None, None)


def _model(plddt: list) -> Model:
    return Model("A" * len(plddt), np.array(plddt, dtype=float),
                 np.array([[3.8 * i, 0.0, 0.0] for i in range(len(plddt))]))


def _sln_job(**changes) -> Job:
    return Job(**{**dict(slug="sln", accession="O00631", display="Sarcolipin", protein=SLN,
                         kept=((1, 31),)), **changes})


# ------------------------------------------------------------ the entry


def test_the_canonical_entry_is_taken_from_among_the_isoforms():
    answer = _answer("Q9H6X2")
    assert [e["entryId"] for e in answer][-1] == "AF-Q9H6X2-F1"   # isoforms come first here
    entry = alphafold.entry_of("Q9H6X2", answer, 564)
    assert (entry.entry_id, entry.version) == ("AF-Q9H6X2-F1", 6)
    assert entry.model_url == "https://alphafold.ebi.ac.uk/files/AF-Q9H6X2-F1-model_v6.pdb"


def test_an_answer_of_isoforms_alone_is_no_model():
    isoforms = [e for e in _answer("Q9H6X2") if e["uniprotAccession"] != "Q9H6X2"]
    with pytest.raises(Refused, match="holds no model of UniProt Q9H6X2"):
        alphafold.entry_of("Q9H6X2", isoforms, 564)
    with pytest.raises(Refused, match="holds no model of UniProt Q9H6X2"):
        alphafold.entry_of("Q9H6X2", [], 564)


def test_a_protein_alphafold_db_does_not_hold_is_refused_and_told_why():
    url = alphafold.API.format(accession="O95714")
    get = _serving({url: _not_found(url)})
    with pytest.raises(Refused) as long:
        alphafold.fetch_entry("O95714", 4834, get)
    assert str(long.value) == ("AlphaFold DB has no model of proteins over 2,700 residues; "
                               "this one has 4,834.")
    with pytest.raises(Refused) as short:
        alphafold.fetch_entry("O95714", 300, get)
    assert str(short.value) == "AlphaFold DB holds no model of UniProt O95714."


def test_an_api_that_is_down_is_not_a_refusal():
    url = alphafold.API.format(accession="O00631")
    down = urllib.error.HTTPError(url, 503, "Service Unavailable", None, None)
    with pytest.raises(urllib.error.HTTPError):
        alphafold.fetch_entry("O00631", 31, _serving({url: down}))


def test_the_model_is_read_residue_by_residue():
    model = alphafold.read_model(SLN_MODEL.read_text())
    assert model.letters == SLN
    assert model.ca.shape == (31, 3)
    # The API's own mean over the whole chain, to the tenth the page shows it
    # to: the file rounds each residue's pLDDT to two decimals, and the mean
    # of those is 91.57 where the API says 91.56.
    assert abs(float(model.plddt.mean()) - _answer("O00631")[0]["globalMetricValue"]) < 0.05


def test_a_file_that_is_not_one_plain_chain_is_not_read_as_a_model():
    lines = [line for line in SLN_MODEL.read_text().splitlines() if line.startswith("ATOM")]
    gapped = [line for line in lines if int(line[22:26]) != 10]
    with pytest.raises(ValueError, match="residue 10 is numbered 11"):
        alphafold.read_model("\n".join(gapped))
    other = [line[:21] + "B" + line[22:] for line in lines]
    with pytest.raises(ValueError, match="not one plain chain A"):
        alphafold.read_model("\n".join(other))
    with pytest.raises(ValueError, match="holds no residues"):
        alphafold.read_model("HEADER\n")


def test_the_model_is_held_to_the_records_protein():
    model = alphafold.read_model(SLN_MODEL.read_text())
    assert alphafold.check_sequence("O00631", model, SLN, 0) == 0
    one_off = "A" + SLN[1:]
    assert alphafold.check_sequence("O00631", model, one_off, 3) == 1
    with pytest.raises(Refused, match="differs from this record's protein at 1 residues"):
        alphafold.check_sequence("O00631", model, one_off, 0)
    with pytest.raises(Refused) as longer:
        alphafold.check_sequence("O00631", model, SLN + "K", 3)
    assert str(longer.value) == ("AlphaFold's model of O00631 is 31 residues long; "
                                 "this record's protein is 32.")


# ------------------------------------------------------------ what is drawn


def test_the_span_is_the_first_kept_region_to_the_last():
    assert alphafold.span_of(((21, 119),), 119) == (21, 119)           # a leader is left out
    assert alphafold.span_of(((25, 54), (90, 110)), 110) == (25, 110)  # chains, and what is between
    assert alphafold.span_of((), 838) == (1, 838)                      # no region named
    with pytest.raises(ValueError):
        alphafold.span_of(((5, 120),), 119)


def test_the_bands_are_alphafolds_own():
    assert [name for name, _ in BANDS] == [
        "plddtVeryHigh", "plddtConfident", "plddtLow", "plddtVeryLow"]
    scores = np.array([95.0, 90.0, 89.99, 70.0, 69.99, 50.0, 49.99, 3.0])
    assert alphafold.band_of(scores).tolist() == [0, 0, 1, 1, 2, 2, 3, 3]
    # Counted this way, a whole model's shares are the ones the API states.
    model = alphafold.read_model(SLN_MODEL.read_text())
    stated = _answer("O00631")[0]
    assert [round(s, 3) for s in alphafold.confidence(model, (1, 31)).shares] == [
        stated["fractionPlddtVeryHigh"], stated["fractionPlddtConfident"],
        stated["fractionPlddtLow"], stated["fractionPlddtVeryLow"]]


def test_a_model_under_the_gate_is_refused_with_its_number():
    unsure = _model([41.6] * 1863)
    said = alphafold.confidence(unsure, (1, 1863))
    with pytest.raises(Refused) as refused:
        alphafold.gate(said, (1, 1863))
    assert str(refused.value) == (
        "AlphaFold's model is not confident here: mean pLDDT 41.6 over residues "
        "1–1,863, under the 50 needed to draw it.")
    alphafold.gate(alphafold.confidence(_model([50.0] * 10), (1, 10)), (1, 10))


def test_the_gate_reads_the_span_and_not_the_leader():
    # A signal peptide the model is unsure of does not refuse a chain it is sure of.
    model = _model([20.0] * 20 + [90.0] * 10)
    assert alphafold.confidence(model, (1, 30)).mean < alphafold.GATE
    said = alphafold.confidence(model, (21, 30))
    assert said.mean == 90.0 and said.shares == (1.0, 0.0, 0.0, 0.0)
    alphafold.gate(said, (21, 30))


def test_a_bridge_is_drawn_only_where_the_models_sulfurs_meet():
    at = {
        ("A", 45, "SG"): np.array([0.0, 0.0, 0.0]), ("A", 100, "SG"): np.array([2.02, 0.0, 0.0]),
        ("A", 60, "SG"): np.array([9.0, 0.0, 0.0]), ("A", 70, "SG"): np.array([9.0, 5.3, 0.0]),
        ("A", 5, "SG"): np.array([1.0, 1.0, 1.0]), ("A", 50, "SG"): np.array([1.0, 1.0, 3.0]),
    }
    drawn, dropped = alphafold.bridges_of(
        ((100, 45), (60, 70), (5, 50), (80, 90)), at, (21, 119))
    assert drawn == [(45, 100)]
    assert dropped == [
        {"pair": [5, 50], "reason": "outside the span drawn"},
        {"pair": [60, 70], "reason": "sulfurs apart in the model", "distance": 5.3},
        {"pair": [80, 90], "reason": "no sulfur in the model"},
    ]


def test_the_bakes_copy_names_its_bridges_as_an_entry_does(tmp_path):
    path = alphafold.working_copy(SLN_MODEL.read_text(), (5, 20), [(7, 15), (9, 18)],
                                  tmp_path / "span.pdb")
    structure = Structure(pdb="AF-O00631-F1", chains=(Chain("chainA", "A"),), bonds=True)
    # Read back by the bake's own reader of an entry's SSBOND records.
    assert ssbonds(path, structure) == [("A", 7, "A", 15), ("A", 9, "A", 18)]
    numbers = {int(line[22:26]) for line in path.read_text().splitlines()
               if line.startswith("ATOM")}
    assert numbers == set(range(5, 21))
    assert len(ca_atoms(structure, path)) == 16


def test_sampling_falls_as_the_span_grows():
    assert alphafold.sampling_for(31) == 8
    assert alphafold.sampling_for(99) == 8
    assert alphafold.sampling_for(838) == 1
    assert alphafold.sampling_for(2700) == 1
    chosen = [alphafold.sampling_for(n) for n in range(1, 3001)]
    assert all(a >= b for a, b in zip(chosen, chosen[1:]))
    assert set(chosen) <= set(range(1, 9))


def test_per_cents_add_to_a_hundred():
    assert alphafold.percentages((0.774, 0.161, 0.065, 0.0)) == [77, 16, 7, 0]
    assert alphafold.percentages((1 / 3, 1 / 3, 1 / 3, 0.0)) == [34, 33, 33, 0]
    assert alphafold.percentages((1.0, 0.0, 0.0, 0.0)) == [100, 0, 0, 0]
    for shares in ((0.2955, 0.4195, 0.0865, 0.1985), (0.005, 0.005, 0.495, 0.495)):
        assert sum(alphafold.percentages(shares)) == 100


def test_the_fold_pages_words_are_the_models_own_numbers():
    said = Confidence(96.97, (0.96, 0.03, 0.01, 0.0))
    job = Job(slug="b2m", accession="P61769", display="Beta-2-microglobulin",
              protein="M" * 119, kept=((21, 119),), disulfides=((45, 100),))
    entry = Entry("P61769", "AF-P61769-F1", 6, "https://example.invalid/model.pdb")
    chrome, chains = alphafold.describe(
        job, entry, (21, 119), said, [(45, 100)],
        ["plddtVeryHigh", "plddtConfident", "plddtLow", "bonds"])
    # The seven keys the twenty's rows carry, and no other.
    assert list(chrome) == ["pdb", "modelled", "label", "count", "unit", "sentence", "semantics"]
    assert chrome["pdb"] == "AF-P61769-F1"
    assert chrome["modelled"] == [21, 119]
    assert (chrome["label"], chrome["count"], chrome["unit"]) == (
        "the predicted fold", 99, "residues")
    assert chrome["sentence"] == (
        "AlphaFold prediction, mean pLDDT 97.0. 99% of residues at 70 or over.")
    assert chrome["semantics"] == (
        "AlphaFold's predicted fold of Beta-2-microglobulin, residues 21 to 119, coloured by "
        "confidence: 96% very high, 3% confident, 1% low. "
        "Its disulfide bridge is drawn. Drag to turn it.")
    assert chains == [
        {"node": "plddtVeryHigh", "tint": "plddtVeryHigh"},
        {"node": "plddtConfident", "tint": "plddtConfident"},
        {"node": "plddtLow", "tint": "plddtLow"},
        {"node": "bonds", "tint": "cysteine"},
    ]


def test_a_model_of_the_whole_protein_names_no_span_and_groups_its_thousands():
    job = Job(slug="mtor", accession="P42345", display="mTOR", protein="M" * 2549,
              kept=((1, 2549),))
    entry = Entry("P42345", "AF-P42345-F1", 6, "https://example.invalid/model.pdb")
    chrome, chains = alphafold.describe(job, entry, (1, 2549),
                                        Confidence(78.0, (0.4, 0.4, 0.1, 0.1)), [], ["plddtLow"])
    assert chrome["modelled"] is None
    assert chrome["count"] == 2549
    assert "residues 1 to 2,549" in chrome["semantics"]
    assert "bridge" not in chrome["semantics"]
    assert chains == [{"node": "plddtLow", "tint": "plddtLow"}]


def test_only_the_bands_a_model_has_are_spoken():
    job = Job(slug="x", accession="P00000", display="X", protein="M" * 400, kept=((1, 400),))
    entry = Entry("P00000", "AF-P00000-F1", 6, "https://example.invalid/model.pdb")
    # One residue in four hundred is a quarter of a per cent, and still there.
    chrome, _ = alphafold.describe(job, entry, (1, 400),
                                   Confidence(88.0, (0.0, 0.9975, 0.0025, 0.0)), [], [])
    assert "coloured by confidence: 100% confident, under 1% low. Drag" in chrome["semantics"]


def test_a_cut_precursor_is_said_to_be_drawn_as_one_chain():
    job = Job(slug="x", accession="P00000", display="X", protein="M" * 110,
              kept=((25, 54), (90, 110)))
    entry = Entry("P00000", "AF-P00000-F1", 6, "https://example.invalid/model.pdb")
    chrome, _ = alphafold.describe(job, entry, (25, 110),
                                   Confidence(62.14, (0.1, 0.4, 0.3, 0.2)), [], ["plddtLow"])
    assert chrome["sentence"] == (
        "AlphaFold's precursor, mean pLDDT 62.1: every chain and what lies between them.")
    assert chrome["modelled"] == [25, 110]


def test_no_sentence_outgrows_its_two_lines():
    # R4.5: about 80 characters on a phone, at the longest numbers there are.
    entry = Entry("P00000", "AF-P00000-F1", 6, "https://example.invalid/model.pdb")
    for kept in (((1, 2700),), ((1, 1000), (1200, 2700))):
        job = Job(slug="x", accession="P00000", display="X", protein="M" * 2700, kept=kept)
        chrome, _ = alphafold.describe(job, entry, (1, 2700),
                                       Confidence(100.0, (1.0, 0.0, 0.0, 0.0)), [], [])
        assert len(chrome["sentence"]) <= 80, chrome["sentence"]


def test_a_scene_is_compiled_only_by_the_importer_the_stored_scenes_were(tmp_path):
    def lock(version: str) -> Path:
        (tmp_path / "pubspec.lock").write_text(
            'packages:\n  flutter_bloc:\n    dependency: "direct main"\n    version: "9.1.1"\n'
            '  flutter_scene:\n    dependency: "direct main"\n    description:\n'
            '      name: flutter_scene\n      url: "https://pub.dev"\n    source: hosted\n'
            f'    version: "{version}"\n  flutter_test:\n    version: "0.0.0"\n')
        return tmp_path

    assert alphafold.importer_version(lock(alphafold.FLUTTER_SCENE)) == alphafold.FLUTTER_SCENE
    assert alphafold.importer_unready(tmp_path) is None
    assert "pins flutter_scene 0.24.0" in alphafold.importer_unready(lock("0.24.0"))
    assert "names no flutter_scene" in alphafold.importer_unready(tmp_path / "nowhere")
    with pytest.raises(RuntimeError, match="pins flutter_scene 0.24.0"):
        alphafold.compile_scene(tmp_path / "x.glb", tmp_path / "x.fsceneb", tmp_path, "dart")


# ------------------------------------------------------------ the mesh


def test_each_point_of_a_ribbon_goes_to_the_nearer_of_its_two_residues():
    pytest.importorskip("scipy")
    ca = np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0]])
    points = np.array([[-1.0, 0.5, 0.0], [1.8, 0.4, 0.0], [2.0, 0.4, 0.0],
                       [5.6, -0.3, 0.2], [5.8, -0.3, 0.2], [9.0, 0.0, 0.0]])
    assert alphafold.residue_of(points, ca).tolist() == [0, 0, 1, 1, 2, 2]
    assert alphafold.residue_of(points, ca[:1]).tolist() == [0] * 6


def test_a_strands_ribbon_is_shared_out_by_where_it_lies_along_the_chain():
    pytest.importorskip("scipy")
    # A flat arrow is smoothed past the atoms it is drawn for: the CA atoms
    # zigzag a good angstrom to either side of a ribbon that runs straight.
    ca = np.array([[3.3 * i, 1.0 if i % 2 else -1.0, 0.0] for i in range(6)])
    along = np.linspace(0.0, 16.5, 34)
    ribbon = np.stack([along, np.zeros_like(along), np.zeros_like(along)], axis=1)
    owner = alphafold.residue_of(ribbon, ca)
    assert (np.diff(owner) >= 0).all()                     # in order along the chain
    assert np.bincount(owner, minlength=6).min() >= 3      # and every residue has its share


def test_the_bands_are_the_cartoon_triangle_for_triangle():
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    ca = np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [7.6, 0.0, 0.0], [11.4, 0.0, 0.0]])
    origin = np.array([5.0, -2.0, 1.0])
    tube = trimesh.creation.cylinder(radius=0.4, segment=[ca[0], ca[-1]], sections=8)
    tube = tube.subdivide().subdivide()
    tube.apply_translation(-origin)                        # the export frame
    normals = np.asarray(tube.vertex_normals).copy()
    band = np.array([0, 2, 2, 0])                          # no residue in bands 1 and 3
    meshes = alphafold.split_by_band(tube, origin, ca, band)

    assert list(meshes) == ["plddtVeryHigh", "plddtLow"]
    assert sum(len(m.faces) for m in meshes.values()) == len(tube.faces)
    low = meshes["plddtLow"]
    centres = np.asarray(low.vertices)[np.asarray(low.faces)].mean(axis=1) + origin
    assert (centres[:, 0] > 1.9 - 1e-9).all() and (centres[:, 0] < 9.5 + 1e-9).all()
    # Vertices and PyMOL's normals go with their triangles, as they were.
    where = {tuple(np.round(v, 9)): n for v, n in zip(np.asarray(tube.vertices), normals)}
    for mesh in meshes.values():
        for vertex, normal in zip(np.asarray(mesh.vertices), np.asarray(mesh.vertex_normals)):
            assert np.allclose(where[tuple(np.round(vertex, 9))], normal)


def test_a_turned_or_shifted_frame_is_caught():
    pytest.importorskip("scipy")
    ca = np.array([[3.8 * i, 0.0, 0.0] for i in range(12)])
    origin = np.array([4.0, 1.0, -2.0])
    ribbon = np.repeat(ca, 4, axis=0) + np.tile(
        [[0.0, 0.2, 0.0], [0.0, -0.2, 0.0], [0.0, 0.0, 0.2], [0.0, 0.0, -0.2]], (12, 1))
    found = alphafold.frame_audit(ribbon - origin, origin, ca)
    assert found == {"ca_mean": 0.2, "ca_max": 0.2, "ca_on_ribbon": 1.0}
    with pytest.raises(ValueError, match="export frame is off"):
        alphafold.frame_audit(ribbon - origin, origin + np.array([0.0, 10.0, 0.0]), ca)
    with pytest.raises(ValueError, match="export frame is off"):
        alphafold.frame_audit((ribbon - origin)[:, [1, 0, 2]], origin, ca)
    # A chain drawn as a tube lies on the tube's axis, its radius from the surface.
    tube = np.repeat(ca, 4, axis=0) + np.tile(
        [[0.0, 0.6, 0.0], [0.0, -0.6, 0.0], [0.0, 0.0, 0.6], [0.0, 0.0, -0.6]], (12, 1))
    with pytest.raises(ValueError, match="0% of them on it"):
        alphafold.frame_audit(tube - origin, origin, ca)
    assert alphafold.frame_audit(tube - origin, origin, ca, on=0.6)["ca_on_ribbon"] == 1.0


# ------------------------------------------------------------ the whole bake


def _bake_tools():
    """PyMOL, Dart and the app's checkout, or why the whole bake cannot run here."""
    pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    from pipeline.structure import bake

    if not Path(bake.PYMOL).exists():
        pytest.skip("PyMOL is not installed")
    dart = shutil.which("dart") or str(Path.home() / "flutter" / "bin" / "dart")
    if not Path(dart).exists():
        pytest.skip("no Dart SDK to run the scene importer with")
    unready = alphafold.importer_unready(CLIENT)
    if unready:
        pytest.skip(unready)
    return dart


def _sln_pages() -> dict:
    return {
        alphafold.API.format(accession="O00631"): (FIXTURES / "O00631.json").read_bytes(),
        SLN_ENTRY.model_url: SLN_MODEL.read_bytes(),
    }


def test_sarcolipin_is_baked_whole(tmp_path):
    dart = _bake_tools()
    from pipeline.structure.glb import read_glb

    built = alphafold.build(_sln_job(), tmp_path, CLIENT, dart, _serving(_sln_pages()))

    meshes = read_glb(built.glb)
    assert list(meshes) == ["plddtVeryHigh", "plddtConfident", "plddtLow"]
    assert [c["node"] for c in built.chains] == list(meshes)
    everything = np.vstack([m.positions for m in meshes.values()])
    assert abs(np.ptp(everything, axis=0).max() - 1.0) < 1e-6      # cut to the twenty's measure
    assert 0 < len(built.scene) <= 900_000
    assert built.chrome["count"] == 31 and built.chrome["modelled"] is None
    said = built.provenance
    assert (said["entry"], said["model_version"], said["licence"]) == (
        "AF-O00631-F1", 6, "CC BY 4.0")
    assert said["sampling"] == 8 and said["nodes"] == list(meshes)
    assert said["mean_plddt"] == 91.57
    assert said["frame"]["ca_max"] <= alphafold.FRAME_NEAR
    assert said["importer"] == {"package": "flutter_scene", "version": alphafold.FLUTTER_SCENE}


def test_a_model_no_sampling_fits_in_the_budget_is_refused(tmp_path):
    dart = _bake_tools()
    with pytest.raises(Refused) as refused:
        alphafold.build(_sln_job(), tmp_path, CLIENT, dart, _serving(_sln_pages()), budget=1000)
    assert str(refused.value) == "AlphaFold's model of 31 residues is too large to draw."


def test_a_refusal_comes_before_any_bake(tmp_path):
    # Nothing here needs PyMOL: a protein that is refused costs two requests.
    pages = _sln_pages()
    with pytest.raises(Refused, match="31 residues long; this record's protein is 32"):
        alphafold.build(_sln_job(protein=SLN + "K", kept=((1, 32),)), tmp_path, CLIENT,
                        "dart", _serving(pages))
    assert not list(tmp_path.iterdir())
