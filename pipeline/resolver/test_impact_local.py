"""The AVI bake's driver, offline: the bake called unchanged, its verdicts told
from its failures, and nothing written where the twenty's are kept.

The AlphaGenome client is not in this environment, so where the driver imports
it a stand-in module is put in its place; `bake_impact` itself imports neither
the client nor numpy until it pulls, and nothing here pulls.
"""

from __future__ import annotations

import os
import pickle
import subprocess
import sys
import types
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.impact import bake_impact  # noqa: E402
from pipeline.resolver import impact_local  # noqa: E402
from pipeline.resolver.test_local_worker import INS  # noqa: E402

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="the worker is a Mac's: launchd, flock, POSIX signals")

# A record with three exons and two introns of 180 and 788 bases, and one of a
# single exon: what `intron_spans` and `check_biology` read of a record.
THREE_EXONS = {"location": {"start": 1, "end": 1200, "strand": 1},
               "exons": [{"start": 1, "end": 40}, {"start": 221, "end": 400},
                         {"start": 1189, "end": 1200}]}
ONE_EXON = {"location": {"start": 1, "end": 500, "strand": 1},
            "exons": [{"start": 1, "end": 500}]}


@pytest.fixture
def unguarded(monkeypatch):
    """`bake_impact` as it was imported, put back after each test however the
    driver wrapped or pointed it."""
    for name in ("pull", "gencode", "check_biology", "bake", "DATA", "GENCODE_DIR",
                 "CHECKPOINT_DIR", "DEFAULT_SKILL"):
        monkeypatch.setattr(bake_impact, name, getattr(bake_impact, name))
    monkeypatch.setattr(bake_impact, "_guarded", False, raising=False)
    return bake_impact


@pytest.fixture
def client(monkeypatch):
    """The AlphaGenome client as the driver imports it: `atlas.create(key)`
    hands back what it was given, and a query answers as `answer` says."""
    said = types.SimpleNamespace(answer=lambda interval, scorers: {"AVI_SCORE": "scores"},
                                 asked=[])

    class Client:
        def __init__(self, key):
            self.key = key

        def query_interval(self, interval, requested_scorers):
            said.asked.append((interval, requested_scorers))
            return said.answer(interval, requested_scorers)

    alphagenome = types.ModuleType("alphagenome")
    atlas = types.ModuleType("alphagenome.atlas")
    atlas.atlas = types.SimpleNamespace(create=Client)
    data = types.ModuleType("alphagenome.data")
    data.genome = types.SimpleNamespace(Interval=lambda **where: where)
    for name, module in (("alphagenome", alphagenome), ("alphagenome.atlas", atlas),
                         ("alphagenome.data", data)):
        monkeypatch.setitem(sys.modules, name, module)
    return said


def _hand_over(tmp_path: Path) -> list:
    tmp_path.mkdir(parents=True, exist_ok=True)
    target_file, record_file = tmp_path / "target.pkl", tmp_path / "record.json"
    target_file.write_bytes(pickle.dumps(INS))
    record_file.write_bytes(b'{"gene": "INS"}\n')
    return [str(target_file), str(record_file), str(tmp_path / "out"),
            "--state", str(tmp_path / "state")]


def test_the_bake_reads_and_writes_in_the_workers_folders_and_is_called_unchanged(
        tmp_path, unguarded):
    asked = []

    def bake(target, client, skill, resume, refresh_gencode):
        asked.append((target, client, skill, resume, refresh_gencode, bake_impact.DATA,
                      bake_impact.GENCODE_DIR, bake_impact.CHECKPOINT_DIR,
                      (bake_impact.DATA / target.mock_asset).read_bytes()))
        out = bake_impact.DATA / target.impact_asset
        out.parent.mkdir(parents=True)
        out.write_bytes(b"the track")

    unguarded.bake = bake
    data, state = tmp_path / "data", tmp_path / "state"
    assert impact_local.bake(INS, b"the record", data, state, "client") == b"the track"
    assert asked == [(INS, "client", bake_impact.DEFAULT_SKILL, True, False, data,
                      state / "gencode", state / "checkpoints", b"the record")]
    # Never the folders the twenty's GENCODE answers are committed in.
    here = BACKEND / "pipeline" / "impact"
    assert here not in (data, *data.parents) and here not in (state, *state.parents)


def test_a_gate_is_a_verdict_and_takes_its_checkpoint_with_it(tmp_path, unguarded):
    gate = "INS: record has 3 exons, MANE Select has 4"
    unguarded.bake = lambda *said, **more: (_ for _ in ()).throw(bake_impact.BakeError(gate))
    checkpoint = tmp_path / "state" / "checkpoints" / "ins.jsonl"
    checkpoint.parent.mkdir(parents=True)
    checkpoint.write_text('{"w": "chr11:1-2", "d": {}}\n')
    with pytest.raises(bake_impact.BakeError, match="MANE Select has 4"):
        impact_local.bake(INS, b"{}", tmp_path / "data", tmp_path / "state", "client")
    assert not checkpoint.exists()


def test_the_atlas_not_answering_is_never_a_verdict(unguarded):
    def pull(client, chromosome, low, high):
        raise bake_impact.BakeError(f"{chromosome}:{low}-{high} failed after 7 attempts: "
                                    f"RESOURCE_EXHAUSTED")

    unguarded.pull = pull
    impact_local.guard(bake_impact)
    with pytest.raises(impact_local.Unreachable, match="RESOURCE_EXHAUSTED") as raised:
        bake_impact.pull("client", "chr11", 1, 2)
    assert not isinstance(raised.value, bake_impact.BakeError)
    # Wrapped once, however often the driver is asked.
    wrapped = bake_impact.pull
    impact_local.guard(bake_impact)
    assert bake_impact.pull is wrapped


def test_gencode_failing_is_tried_again_and_no_single_mane_select_is_a_verdict(unguarded):
    said = {}

    def gencode(gene, skill, refresh=False):
        raise bake_impact.BakeError(said["error"])

    unguarded.gencode = gencode
    impact_local.guard(bake_impact)
    said["error"] = "gtf failed for INS: uv: command not found"
    with pytest.raises(impact_local.Unreachable, match="gtf failed"):
        bake_impact.gencode("INS", Path("/skill"))
    said["error"] = "SHOX: expected exactly one MANE Select transcript, found 2"
    with pytest.raises(bake_impact.BakeError, match="found 2"):
        bake_impact.gencode("SHOX", Path("/skill"))


def test_the_bases_the_biology_alarm_compares_are_counted_as_it_counts_them():
    # 180 and 788 intron bases, less eight at each end of each.
    assert impact_local.interior(bake_impact, THREE_EXONS) == (180 - 16) + (788 - 16)
    assert impact_local.interior(bake_impact, ONE_EXON) == 0
    # An intron of ten has an edge of five at each end, and so no interior.
    short = {**THREE_EXONS, "exons": [{"start": 1, "end": 40}, {"start": 51, "end": 60}]}
    assert impact_local.interior(bake_impact, short) == 0


def test_a_gene_with_no_intron_skips_only_the_biology_alarm(unguarded):
    asked = []
    unguarded.check_biology = lambda record, positions, gene: asked.append(gene) or "checked"
    impact_local.guard(bake_impact)

    # A gene with introns is checked as the bake checks it.
    assert bake_impact.check_biology(THREE_EXONS, {1: [9.0, 9.0, 9.0]}, "INS") == "checked"
    # One with none has no interior to compare its exons with: its exons'
    # median is kept, and the rest said to be none.
    positions = {local: [float(local % 7), 1.0, 2.0] for local in range(1, 501)}
    assert bake_impact.check_biology(ONE_EXON, positions, "ADRB2") == {
        "exon_median_phred": 3.0, "splice_median_phred": None, "intron_median_phred": None}
    # Nothing scored is still the bake's to say.
    assert bake_impact.check_biology(ONE_EXON, {}, "SRY") == "checked"
    assert asked == ["INS", "SRY"]


def test_the_driver_writes_the_track_and_exits_0(tmp_path, monkeypatch, unguarded, client):
    handed = []

    def bake(target, record, data, state, made):
        handed.append((target, record, data, state, made.key))
        return b"the track"

    monkeypatch.setattr(impact_local, "bake", bake)
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "the-atlas-key")
    monkeypatch.setenv("HELIXPEEK_DATA", "untouched")
    assert impact_local.main(_hand_over(tmp_path)) == 0
    assert (tmp_path / "out" / "impact.json").read_bytes() == b"the track"
    assert handed == [(INS, b'{"gene": "INS"}\n', tmp_path / "out" / "data",
                       tmp_path / "state", "the-atlas-key")]
    assert os.environ["HELIXPEEK_DATA"] == str(tmp_path / "out" / "data")


def test_the_driver_writes_the_gates_words_and_exits_3(tmp_path, monkeypatch, unguarded, client):
    gate = "INS: only 41.2% of drawn bases scored (842 without a full set of three)"
    monkeypatch.setattr(impact_local, "bake", lambda *said: (_ for _ in ()).throw(
        bake_impact.BakeError(gate)))
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "the-atlas-key")
    assert impact_local.main(_hand_over(tmp_path)) == impact_local.REFUSED == 3
    assert (tmp_path / "out" / "refusal.txt").read_text() == gate
    assert not (tmp_path / "out" / "impact.json").exists()

    # The bake breaking is no verdict: it is raised, and Python exits 1.
    monkeypatch.setattr(impact_local, "bake", lambda *said: (_ for _ in ()).throw(
        impact_local.Unreachable("chr11:1-2 failed after 7 attempts: UNAVAILABLE")))
    with pytest.raises(impact_local.Unreachable):
        impact_local.main(_hand_over(tmp_path / "again"))


def test_the_driver_needs_its_key_and_its_four_arguments(tmp_path, monkeypatch, client, capsys):
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    assert impact_local.main(_hand_over(tmp_path)) == 1
    assert "ALPHAGENOME_API_KEY is not set" in capsys.readouterr().err
    assert impact_local.main([str(tmp_path / "target.pkl")]) == 2


def test_the_driver_imports_nothing_heavy_until_it_is_asked():
    # The worker imports it for its exit status, in an environment with no
    # AlphaGenome client; and `pipeline.paths` must wait for HELIXPEEK_DATA.
    code = ("import sys; import pipeline.resolver.impact_local as m; "
            "print(m.REFUSED, sorted(n for n in sys.modules if n.split('.')[0] in "
            "('alphagenome', 'numpy') or n in ('pipeline.paths', 'pipeline.impact')))")
    done = subprocess.run([sys.executable, "-c", code], cwd=BACKEND, capture_output=True,
                          text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "3 []"


@pytest.fixture
def skill(tmp_path, unguarded):
    folder = tmp_path / "alphagenome-variant-impact-score"
    (folder / "scripts").mkdir(parents=True)
    (folder / "scripts" / "alphagenome_atlas_avi.py").write_text("")
    unguarded.DEFAULT_SKILL = folder
    return folder


def test_the_check_asks_the_atlas_once_and_says_what_is_missing(
        skill, monkeypatch, client, capsys):
    ran = []
    looked_up = subprocess.CompletedProcess([], 0, "usage: gtf", "")
    monkeypatch.setattr(impact_local.subprocess, "run",
                        lambda command, **how: ran.append((command, how["cwd"])) or looked_up)

    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    assert impact_local.check() == 1
    assert capsys.readouterr().err == "ALPHAGENOME_API_KEY is not set in the worker's .env.\n"

    monkeypatch.setenv("ALPHAGENOME_API_KEY", "the-atlas-key")
    assert impact_local.check() == 0
    assert ran == [(["uv", "run", str(skill / "scripts" / "alphagenome_atlas_avi.py"), "gtf",
                     "--help"], skill)]
    # One base, one request.
    assert client.asked == [({"chromosome": "chr18", "start": 31_591_766, "end": 31_591_767},
                             ["AVI_SCORE"])]

    client.answer = lambda interval, scorers: (_ for _ in ()).throw(
        ValueError("API key not valid. Please pass a valid API key."))
    assert impact_local.check() == 1
    assert "The Atlas did not answer: API key not valid" in capsys.readouterr().err
    assert "the-atlas-key" not in capsys.readouterr().err

    monkeypatch.setattr(impact_local.subprocess, "run", lambda command, **how: (
        subprocess.CompletedProcess(command, 2, "", "error: No interpreter found for 3.14")))
    assert impact_local.check() == 1
    assert "GENCODE's lookup exited with status 2: error: No interpreter" in \
        capsys.readouterr().err

    (skill / "scripts" / "alphagenome_atlas_avi.py").unlink()
    assert impact_local.check() == 1
    assert "There is no AlphaGenome skill at" in capsys.readouterr().err
