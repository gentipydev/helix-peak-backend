"""Audit an export frame: every CA within 0.25 A of the exported ribbon.

A rotation would show up as tens of angstroms. With no arguments it audits
insulin: run it after `bake.py --target insulin`, which leaves the .obj files
in output/insulin/. `--pdb`, `--out`, `--origin` and `--chain` audit another
export against the file it was made from, as `assemblies/` does for its pair.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.structure.pdb import atoms  # noqa: E402  the bake's own reading

ORIGIN = np.array([-18.0825, -1.3525, -9.326])
PDB = 'structures/3I40.pdb'
OUT = 'output/insulin'
CHAINS = (('A', 'chainA.obj'), ('B', 'chainB.obj'))

def ca_atoms(chain, pdb=PDB):
    return np.array([xyz for (c, _, name), xyz in atoms(Path(pdb)).items()
                     if c == chain and name == 'CA'])

def audit(pdb=PDB, out=OUT, origin=ORIGIN, chains=CHAINS):
    """Each chain's CA atoms against its exported ribbon; prints and returns
    (chain, min, mean, max) in angstroms."""
    import trimesh  # the structure bake's environment, as the export is

    found = []
    for chain, obj in chains:
        m = trimesh.load(f'{out}/{obj}', process=False)
        v = np.asarray(m.vertices) + origin          # export frame -> PDB frame
        ca = ca_atoms(chain, pdb)
        d = np.linalg.norm(ca[:, None, :] - v[None, :, :], axis=2).min(axis=1)
        print(f'chain {chain}: {len(ca)} CA vs {len(v)} verts | '
              f'min {d.min():.2f} A  mean {d.mean():.2f} A  max {d.max():.2f} A')
        found.append((chain, float(d.min()), float(d.mean()), float(d.max())))
    return found

def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--pdb', default=PDB)
    parser.add_argument('--out', default=OUT)
    parser.add_argument('--origin', type=float, nargs=3, default=list(ORIGIN))
    parser.add_argument('--chain', action='append', metavar='CHAIN=OBJ',
                        help='a PDB chain and the .obj its ribbon was exported to')
    args = parser.parse_args(argv)
    chains = (tuple(tuple(c.split('=', 1)) for c in args.chain)
              if args.chain else CHAINS)
    return args.pdb, args.out, np.array(args.origin), chains

if __name__ == '__main__':
    audit(*arguments())
