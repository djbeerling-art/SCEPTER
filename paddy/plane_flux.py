#!/usr/bin/env python3
# Basalt-minus-control fluxes through depth planes from the model's per-cell flux terms (flx/zint_*.txt, v21+).
# Usage: python3 paddy/plane_flux.py --ctrl=<run dir>[,<run dir>...] --trt=<run dir>[,...] [--z=0.2,0.24,bottom]
# The transport terms (adv + dif, incl. dispersion) of the cells above a plane sum to (flux out through the plane) minus
# (inflow at the surface); the surface inflow (irrigation water) is the same in both plots, so the treatment-minus-control
# difference of that sum is the extra flux through the plane, exact and including dispersion and diffusion.
# PDA (potentially durable alkalinity) = ALK - 2 Fe(II) - 2 S(-II) - 2 NH4+ ; sulfide is counted at 2 eq per S (upper
# bound: HS- already in ALK would count 1) and NH4+ (not in the zint files) is omitted here; see alk_flux.py.
import sys, os
import numpy as np, pandas as pd
opt = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--'))
ctrl, trt = opt['ctrl'].split(','), opt['trt'].split(',')
planes = opt.get('z', '0.2,0.24,bottom').split(',')
SP = [('ALK', 'co2sp', 'ALK'), ('Fe2', 'aq', 'fe2'), ('S_T', 'aq', 'hs'), ('Ca', 'aq', 'ca'), ('Mg', 'aq', 'mg')]
def out_above(runs, kind, sp, zp):
    tot = 0.
    for r in runs:
        d = pd.read_csv(os.path.join(r, 'flx', f'zint_{kind}-{sp}.txt'), sep=r'\s+')
        m = d['z'] < zp if zp is not None else np.ones(len(d), bool)
        tot += float((d['adv'] + d['dif'])[m].sum())
    return tot
rows = []
for p in planes:
    zp = None if p == 'bottom' else float(p)
    r = {'plane': p}
    for k, kind, sp in SP:
        r[k] = out_above(trt, kind, sp, zp) - out_above(ctrl, kind, sp, zp)
    r['PDA'] = r['ALK'] - 2*r['Fe2'] - 2*r['S_T']
    r['PDA_tCO2/ha'] = r['PDA']*0.4401
    rows.append(r)
print(pd.DataFrame(rows).to_string(index=False, float_format='%.4f'))
