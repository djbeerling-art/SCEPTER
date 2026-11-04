"""Grain iAs change implied by the modelled root As(III) uptake ratio, using the empirical link across the 15 published
Si treatments (studies.csv): log10(grain ratio) = a + b log10(predicted uptake ratio, Ki 300 uM).
Reports the central estimate, a bootstrap 5-95% range of the fitted line (treatments resampled), and the proportional
bound (b = 1, a = 0: grain follows uptake one-for-one; the fitted b < 1 is partly attenuation from noisy porewater data).
Usage: python3 paddy/arsenic/grain_from_uptake.py as_sweep_summary.csv > grain_from_uptake.csv"""
import sys, os, csv
import numpy as np
Ki = 300.0
rows = [r for r in csv.DictReader(l for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'studies.csv')) if not l.startswith('#'))]
x = np.array([np.log10((float(r['as1'])/float(r['as0']))/((1+float(r['si1_uM'])/Ki)/(1+float(r['si0_uM'])/Ki))) for r in rows])
y = np.array([np.log10(1+float(r['obs_grain_pct'])/100) for r in rows])
b, a = np.polyfit(x, y, 1)
rng = np.random.default_rng(0); boot = []
for _ in range(5000):
    i = rng.integers(0, len(x), len(x))
    if np.ptp(x[i]) > 0: boot.append(np.polyfit(x[i], y[i], 1))
boot = np.array(boot)
print(f'# fit (n={len(x)}): log10 grain = {a:+.3f} + {b:.3f} log10 uptake; r = {np.corrcoef(x, y)[0,1]:.2f}', file=sys.stderr)
w = csv.writer(sys.stdout)
w.writerow(['soil', 'year', 'R_up_Ki300', 'grain_ratio_central', 'grain_ratio_p05', 'grain_ratio_p95', 'grain_ratio_proportional', 'grain_change_pct_central'])
for r in csv.DictReader(open(sys.argv[1])):
    u = np.log10(float(r['R_up_Ki300']))
    g = 10**(boot[:, 1] + boot[:, 0]*u); c = 10**(a + b*u)
    w.writerow([r['soil'], r['year'], r['R_up_Ki300'], round(c, 3), round(np.percentile(g, 5), 3), round(np.percentile(g, 95), 3),
                r['R_up_Ki300'], round(100*(c-1))])
