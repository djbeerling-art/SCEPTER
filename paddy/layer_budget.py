#!/usr/bin/env python3
# Layer budgets from the time-integrated per-cell flux terms (flx/zint_*.txt, model v21+): Fe and carbonate alkalinity.
# Usage: python3 paddy/layer_budget.py <run dir> [<run dir> ...] [--layers=0:0.2,0.2:0.24,0.24:9]
# For each layer and pool: storage (tflx), transport (adv, dif), supply (rain), reactions (other terms), residual (res),
# in mol m-2 over the run(s). SCEPTER's sign convention: every term on the left-hand side, tflx + adv + dif + rxn ... = res,
# with tflx = storage change: a positive transport term is a net LOSS from the layer, a negative rain or reaction term a gain.
# Fe pools: total aqueous Fe(II) and Fe(III) (SCEPTER's aqueous totals include the exchangeable fraction; zint_aqex-fe2
# is that part alone, so it is not added) and the Fe solids (green rust counted as 3 Fe).
import sys, os
import numpy as np, pandas as pd
args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--'))
layers = [tuple(float(x) for x in l.split(':')) for l in opt.get('layers', '0:0.2,0.2:0.24,0.24:9').split(',')]
FE = [('aq', 'fe2', 1), ('aq', 'fe3', 1), ('sld', 'amfe3', 1), ('sld', 'gt', 1), ('sld', 'sid', 1),
      ('sld', 'gr', 3), ('sld', 'fes', 1), ('sld', 'fa', 2)]
TERMS = ['tflx', 'adv', 'dif', 'rain']

def load(kind, sp):
    tot = None
    for run in args:
        f = os.path.join(run, 'flx', f'zint_{kind}-{sp}.txt')
        if not os.path.exists(f): return None
        d = pd.read_csv(f, sep=r'\s+')
        tot = d if tot is None else tot.add(d.drop(columns='z'), fill_value=0).assign(z=d['z'])
    return tot

def summ(d, lo, hi, w=1):
    m = (d['z'] > lo) & (d['z'] < hi); s = d[m].drop(columns='z').sum() * w
    out = {t: s.get(t, 0.) for t in TERMS}
    out['rxn'] = s.drop(labels=[c for c in TERMS + ['res'] if c in s.index]).sum(); out['res'] = s.get('res', 0.)
    return out

for lo, hi in layers:
    print(f'\nlayer {lo:g}-{hi:g} m (mol m-2)')
    print(f'{"pool":14s}' + ''.join(f'{t:>10s}' for t in TERMS + ['rxn', 'res']))
    tot = dict.fromkeys(TERMS + ['rxn', 'res'], 0.)
    for kind, sp, w in FE:
        d = load(kind, sp)
        if d is None: continue
        r = summ(d, lo, hi, w)
        if max(abs(v) for v in r.values()) < 1e-6: continue
        print(f'{kind + ":" + sp:14s}' + ''.join(f'{r[t]:10.4f}' for t in TERMS + ['rxn', 'res']))
        for t in tot: tot[t] += r[t]
    print(f'{"Fe total":14s}' + ''.join(f'{tot[t]:10.4f}' for t in TERMS + ['rxn', 'res']))
    print(f'  Fe: storage change {tot["tflx"]:.4f} = net transport in {-(tot["adv"] + tot["dif"]):.4f} '
          f'(adv {-tot["adv"]:.4f}, dif+disp {-tot["dif"]:.4f}) + supply {-tot["rain"]:.4f} - reactions {tot["rxn"]:.4f} (0 if Fe conserved); '
          f'closure error {tot["tflx"] + tot["adv"] + tot["dif"] + tot["rain"] + tot["rxn"]:.2e}')
    d = load('co2sp', 'ALK')
    if d is not None:
        r = summ(d, lo, hi)
        print('ALK (carbonate):' + ''.join(f' {t} {r[t]:.4f}' for t in TERMS + ['rxn', 'res']))
