"""Year-by-year As(III) uptake ratio (basalt/control) from the model's as3up flux, which the pruned chain years keep
(int_flx_aq-as3.txt), plus season Si ratio from int_flx_aq-si where available.
Usage: python3 paddy/arsenic/as_sweep_years.py <SCEPTER_OUT>/paddy [dose] > as_sweep_years.csv"""
import sys, os, glob, re, csv
import numpy as np
root = sys.argv[1]; dose = sys.argv[2] if len(sys.argv) > 2 else '20'
def tot(run, name, col):
    f = os.path.join(run, 'flx', f'int_flx_{name}.txt')
    try:
        h = open(f).readline().split(); a = np.atleast_2d(np.loadtxt(f, skiprows=1))
        if a.size == 0 or a.shape[1] != len(h) or col not in h: return np.nan
        return a[-1, h.index(col)]*a[-1, 0]
    except (OSError, ValueError): return np.nan
w = csv.writer(sys.stdout); w.writerow(['soil', 'fe3_wt', 'gt_wt', 'irralk', 'year', 'uptake_C_molm2', 'uptake_B_molm2', 'R_up_sim'])
for d0 in sorted(glob.glob(os.path.join(root, 'flood_as_f*_d0_y*s'))):
    m = re.search(r'flood_(as_f([\d.]+)_g([\d.]+)_a([\d.]+))_d0_y(\d+)s$', d0)
    if not m: continue
    tag, f, g, a, y = m.groups(); db = os.path.join(root, f'flood_{tag}_d{dose}_y{y}s')
    uc, ub = tot(d0, 'aq-as3', 'as3up'), tot(db, 'aq-as3', 'as3up')
    w.writerow([tag, f, g, a, y, '%.4g' % uc, '%.4g' % ub, '' if not (np.isfinite(uc) and np.isfinite(ub) and uc) else round(ub/uc, 3)])
