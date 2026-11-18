"""Porewater As response to Si (and pH) in the HFO prototype vs observed porewater As ratios in the uptake-test studies.
Baseline pH values are approximate flooded porewater pH (from each paper where reported, else ~6.5)."""
import sys, csv
sys.path.insert(0, '.')
from hfo_as_si import run
cases = [  # study, treatment, si0 mM, si1 mM, pH0, pH1, observed porewater As ratio, basis
    ('Wang2022', 'husk', 0.206, 0.498, 6.3, 6.65, 1.12, 'total As'),
    ('Wang2022', 'Si-ash', 0.206, 0.705, 6.3, 6.65, 0.88, 'total As'),
    ('Li2009', 'Si gel', 0.011, 0.863, 6.5, 6.5, 1.73, 'As(III)'),
    ('Fleck2013', 'Si pooled', 0.155, 1.121, 6.5, 6.5, 1.23, 'total As'),
    ('SeyfferthFendorf2012', 'Richvale M206 Si-gel', 0.240, 1.570, 6.0, 6.0, 3.6, 'As(III)'),
    ('SeyfferthFendorf2012', 'Richvale IR66 Si-gel', 0.218, 1.430, 6.0, 6.0, 3.38, 'As(III)'),
    ('SeyfferthFendorf2012', 'Hawi M206 Si-gel', 0.050, 1.260, 6.0, 6.0, 2.15, 'As(III)'),
    ('He2020', 'slag 2 t/ha', 0.356, 0.543, 5.5, 6.5, 1.39, 'total As'),
    ('He2020', 'slag 4 t/ha', 0.356, 0.845, 5.5, 7.2, 2.53, 'total As'),
]
import math
def ratios(xs):
    res = []
    for st, tr, s0, s1, p0, p1, obs, basis in cases:
        a, b = run(pH=p0, si_mM=s0, xs=xs), run(pH=p1, si_mM=s1, xs=xs)
        res.append(None if a is None or b is None else (b['As3_uM'] + b['As5_uM']) / (a['As3_uM'] + a['As5_uM']))
    return res
if len(sys.argv) > 1 and sys.argv[1] == 'fit':
    print('xs,log_rmse,n_ok,ratios')
    for xs in (0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
        r = ratios(xs); pairs = [(m, c[6]) for m, c in zip(r, cases) if m]
        e = float('nan') if not pairs else math.sqrt(sum((math.log(m) - math.log(o)) ** 2 for m, o in pairs) / len(pairs))
        print(f"{xs},{e:.3f},{len(pairs)},{' '.join(f'{m:.2f}' if m else 'fail' for m in r)}")
    sys.exit()
w = csv.writer(sys.stdout)
w.writerow(['study', 'treatment', 'Si_fold', 'pH_change', 'model_As3_ratio', 'model_totAs_ratio', 'obs_ratio', 'obs_basis'])
for st, tr, s0, s1, p0, p1, obs, basis in cases:
    a, b = run(pH=p0, si_mM=s0), run(pH=p1, si_mM=s1)
    if a is None or b is None:
        w.writerow([st, tr, round(s1 / s0, 1), p1 - p0, 'fail', 'fail', obs, basis]); continue
    r3 = b['As3_uM'] / a['As3_uM']; rt = (b['As3_uM'] + b['As5_uM']) / (a['As3_uM'] + a['As5_uM'])
    w.writerow([st, tr, round(s1 / s0, 1), round(p1 - p0, 2), round(r3, 2), round(rt, 2), obs, basis])
