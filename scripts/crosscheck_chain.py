#!/usr/bin/env python
"""Re-derive each headline result by an independent route and report both.

Round 18's lesson generalised: the useful check is the one that does not
reuse the derivation being checked. A test written from the same
understanding as the code shares its blind spots.

Usage::

    python scripts/crosscheck_chain.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

import numpy as np
from piezo1.core import Structure
from piezo1.config import STRUCTURE_DIR
from piezo1.structure.geometry import measure_dome
from piezo1.structure.superpose import detect_c3_axis, match_protomers, kabsch
from piezo1.structure.protomers import protomer_blocks
from piezo1.physics.anm import ANM
from piezo1.physics.kinetics import GatingModel
from piezo1.analysis.crosscheck import (dome_curvature_by_parabola,
    dome_curvature_by_cap_geometry, gating_overlap_by_distances,
    t50_by_ode_integration, CrossCheck, compare)
from test_geometry import _tm_surface

checks=[]
st = Structure.from_file(STRUCTURE_DIR/'7WLT.cif')
blocks,_ = protomer_blocks(st); surface = _tm_surface(st,'mouse')
axis = detect_c3_axis(blocks)
primary_R = measure_dome(blocks, surface).radius_of_curvature/10.0
checks.append(CrossCheck('dome radius (exact cap geometry)', primary_R,
    dome_curvature_by_cap_geometry(surface, axis), 'nm', 0.10,
    'algebraic sphere fit', 'exact cap relation R = -(h^2+r^2)/2h, per point'))
c = CrossCheck('dome radius (parabola)', primary_R,
    dome_curvature_by_parabola(surface, axis), 'nm', 0.10,
    'algebraic sphere fit', 'parabola h = h0 - r^2/2R (SMALL-SLOPE)')
c.note = ('expected to disagree: the parabola is a shallow-cap approximation '
          'and PIEZO1 sits at 63 deg. On synthetic caps of known R it is 0.6% '
          'low at 8.6 deg and 25.8% low at 63.4 deg.')
checks.append(c)

flat = Structure.from_file(STRUCTURE_DIR/'7WLU.cif')
_c,cr = protomer_blocks(st); _f,fr = protomer_blocks(flat)
common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))
def resample(s):
    out=[]
    for ch in s.chains:
        m=s.mask_ca()&(s.chain==ch)
        if m.sum()<300: continue
        idx={int(r):i for i,r in enumerate(s.res_seq[m])}; xyz=s.xyz[m]
        if all(r in idx for r in common): out.append(np.array([xyz[idx[r]] for r in common],float))
    return out[:3]
cb, fb = resample(st), resample(flat)
fb=[fb[i] for i in match_protomers(cb,fb).order]
anm = ANM.from_trimer(cb, cutoff=15.0, spring='inverse_square').build()
modes = anm.calc_modes(n_modes=30); anm.label_symmetry(modes)
R,tr,c = kabsch(np.vstack(fb), np.vstack(cb))
disp = (((np.vstack(fb)-c)@R.T+tr) - np.vstack(cb)).ravel()
ov = np.abs(np.asarray(modes.overlap(disp),float))
sym = np.array([s=='A' for s in modes.symmetry])
checks.append(CrossCheck('gating overlap, best A mode', float(ov[sym].max()),
    gating_overlap_by_distances(np.vstack(cb), np.vstack(fb), modes), '', 0.20,
    'Kabsch superposition + Cartesian projection',
    'pairwise distance changes; no superposition'))

m = GatingModel()
checks.append(CrossCheck('half-activation tension', m.half_activation(),
    t50_by_ode_integration(m), 'mN/m', 0.05,
    'matrix exponential, peak of the transient',
    'adaptive Runge-Kutta integration of the master equation'))

print('ADVERSARIAL CROSS-CHECK — each result re-derived by an independent route\n')
compare(checks)
print()
for c in checks:
    print(f"  {c.quantity}:")
    print(f"     primary     = {c.primary_route}")
    print(f"     alternative = {c.alternative_route}")


import json
from piezo1.config import DERIVED_DIR
DERIVED_DIR.mkdir(parents=True, exist_ok=True)
(DERIVED_DIR / "crosscheck.json").write_text(json.dumps(
    [{"quantity": c.quantity, "primary": c.primary,
      "alternative": c.alternative, "relative": c.relative,
      "agrees": c.agrees, "primary_route": c.primary_route,
      "alternative_route": c.alternative_route, "note": c.note}
     for c in checks], indent=1))
print(f"\nwritten to {DERIVED_DIR / 'crosscheck.json'}")
