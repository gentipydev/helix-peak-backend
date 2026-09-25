"""Audit insulin's export frame: every CA within 0.25 A of the exported ribbon.

A rotation would show up as tens of angstroms. Run after `bake.py --target
insulin`, which leaves the .obj files in output/insulin/.
"""

import numpy as np, trimesh

ORIGIN = np.array([-18.0825, -1.3525, -9.326])
PDB = 'structures/3I40.pdb'
OUT = 'output/insulin'

def ca_atoms(chain):
    out = []
    for line in open(PDB):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA' and line[21] == chain:
            out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(out)

for chain, obj in (('A', f'{OUT}/chainA.obj'), ('B', f'{OUT}/chainB.obj')):
    m = trimesh.load(obj, process=False)
    v = np.asarray(m.vertices) + ORIGIN          # export frame -> PDB frame
    ca = ca_atoms(chain)
    d = np.linalg.norm(ca[:, None, :] - v[None, :, :], axis=2).min(axis=1)
    print(f'chain {chain}: {len(ca)} CA vs {len(v)} verts | '
          f'min {d.min():.2f} A  mean {d.mean():.2f} A  max {d.max():.2f} A')
