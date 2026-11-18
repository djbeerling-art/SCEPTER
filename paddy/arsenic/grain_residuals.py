"""Structure of grain residuals against the competitive Si-As(III) uptake model (studies.csv).
residual = log10(observed grain ratio) - log10(predicted uptake ratio), uptake ratio = R_As / T, T = (1+Si1/Ki)/(1+Si0/Ki).
Grouped by treatment type, porewater Si, growth-stage window and As basis (porewater species, grain species), plus the
straw-to-grain partitioning signal (log grain ratio - log straw ratio) as the translocation / node-sequestration test.
No global fit: each grouping is descriptive. Usage: python3 paddy/arsenic/grain_residuals.py [Ki uM, default 300]"""
import sys, csv, os, math
import numpy as np
from scipy import stats
Ki = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
rows = [r for r in csv.DictReader(l for l in open(os.path.join(os.path.dirname(__file__), 'studies.csv')) if not l.startswith('#'))]
def ttype(t):
    t = t.lower()
    if 'gel' in t: return 'Si gel'
    if 'slag' in t: return 'slag (basalt analogue)'
    if 'diatom' in t: return 'diatomite'
    return 'Si-rich amendment (husk, ash, char, silicate)'
def stage(w):
    w = w.lower()
    return 'reproductive only' if ('heading' in w or 'maturity' in w) else 'whole season'
out = []
for r in rows:
    s0, s1, a0, a1 = (float(r[k]) for k in ('si0_uM', 'si1_uM', 'as0', 'as1'))
    T = (1 + s1/Ki)/(1 + s0/Ki); up = (a1/a0)/T
    g = 1 + float(r['obs_grain_pct'])/100
    st = 1 + float(r['obs_straw_pct'])/100 if r['obs_straw_pct'] else None
    out.append(dict(study=r['study'], trt=r['treatment'], type=ttype(r['treatment']), stage=stage(r['window']),
                    pw='As(III)' if 'As(III)' in r['as_basis'] else 'total As',
                    grain='iAs' if ('iAs' in r['grain_basis'] or 'As(III)' in r['grain_basis']) else 'total As',
                    lsi=math.log10(s1/s0), si1=s1, lup=math.log10(up), lg=math.log10(g), res=math.log10(g) - math.log10(up),
                    lst=math.log10(st) if st else None))
w = csv.writer(sys.stdout)
w.writerow(['study', 'treatment', 'type', 'stage', 'porewater_As', 'grain_As', 'Si1_uM', 'log10_RSi', 'log10_uptake_pred',
            'log10_grain_obs', 'residual', 'log10_straw_obs', 'grain_minus_straw'])
for o in out:
    w.writerow([o['study'], o['trt'], o['type'], o['stage'], o['pw'], o['grain'], round(o['si1']), round(o['lsi'], 3),
                round(o['lup'], 3), round(o['lg'], 3), round(o['res'], 3), '' if o['lst'] is None else round(o['lst'], 3),
                '' if o['lst'] is None else round(o['lg'] - o['lst'], 3)])
res = np.array([o['res'] for o in out])
print(f'\n# Ki = {Ki} uM; n = {len(out)}; residual mean {res.mean():+.3f}, sd {res.std(ddof=1):.3f} (log10)', file=sys.stderr)
for key in ('type', 'stage', 'pw', 'grain'):
    print(f'# by {key}:', file=sys.stderr)
    for k in sorted(set(o[key] for o in out)):
        v = np.array([o['res'] for o in out if o[key] == k])
        print(f'#   {k:48s} n={len(v):2d} mean {v.mean():+.3f}' + (f' sd {v.std(ddof=1):.3f}' if len(v) > 1 else ''), file=sys.stderr)
for key, lab in (('lsi', 'log10 R_Si'), ('si1', 'treatment porewater Si')):
    rho, p = stats.spearmanr([o[key] for o in out], res)
    print(f'# Spearman residual vs {lab}: rho {rho:+.2f}, p {p:.3f}', file=sys.stderr)
pairs = [(o['lg'] - o['lst'], o['type'], o['trt']) for o in out if o['lst'] is not None]
print(f'# grain minus straw (log10; >0 = grain fell less than straw): n={len(pairs)}', file=sys.stderr)
for d, t, n in pairs: print(f'#   {n:34s} {t:46s} {d:+.3f}', file=sys.stderr)
