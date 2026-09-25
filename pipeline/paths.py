"""Where the pipeline reads and writes, apart from its own code.

These scripts were `helix-peek/tool/` until Phase 3 of HANDOFF-ONDEMAND.md, and
every path hung off the client's root, because what they baked was the app's
own assets. The app ships no protein data any more: storage holds the tracks,
and a bake writes into a data directory of its own before `upload_tracks.py`
sends it there. `fetch_tracks.py` fills the same directory back from storage
for the tools that read files -- the ClinVar bake reads the record and the AVI
map, `check_assets.py` reads everything.

The directory keeps the `assets/...` layout the client had, so the relative
paths `targets.py` hands out (`assets/mock/gene_ins.json`) resolve under it
unchanged, and a path in a bake's log reads the way it always has.

Each location can be moved with an environment variable. The defaults are the
layout on the dev Mac, where the two repos are siblings.
"""

from __future__ import annotations

import os
from pathlib import Path

PIPELINE = Path(__file__).resolve().parent
BACKEND = PIPELINE.parent
WORKSPACE = BACKEND.parent

# Baked and fetched tracks, in the client's old `assets/...` layout. Gitignored.
DATA = Path(os.environ.get("HELIXPEEK_DATA") or PIPELINE / "data")

# The client checkout, read only to prove its asset files equal storage before
# they are moved out of the bundle.
CLIENT = Path(os.environ.get("HELIXPEEK_CLIENT") or WORKSPACE / "helix-peek")

# The AlphaGenome skills: GENCODE lookups for the impact bake, and the Atlas
# client for the explanations bake.
SKILLS = Path(os.environ.get("ALPHAGENOME_SKILLS") or WORKSPACE / ".claude" / "skills")

# The twenty curated proteins' hand-written fields: display name, summary, chain
# name, facts, chain tints, the fold page's prose, and which families each has.
# They lived in `protein_catalog.dart` until Phase 3; this file is where they are
# edited now.
CURATED = PIPELINE / "curated" / "catalog.json"
