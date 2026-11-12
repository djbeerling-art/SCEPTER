"""Threshold test (David, 2 Nov 2026): does the treatment raise porewater As by more than the Si-inhibition threshold
T = (1 + Si_trt/Ki)/(1 + Si_ctrl/Ki)?  Predicted As(III) uptake ratio = R_As / T. Ki = 300 uM (siup Km) and 150 uM."""
import csv, math
rows = [r for r in csv.DictReader(l for l in open('studies.csv') if not l.startswith('#'))]
out = []
for r in rows:
    s0, s1, a0, a1 = (float(r[k]) for k in ('si0_uM','si1_uM','as0','as1'))
    ra = a1/a0
    rec = dict(study=r['study'], treatment=r['treatment'], si_ratio=round(s1/s0,2), as_ratio=round(ra,2), as_basis=r['as_basis'])
    for ki in (300, 150):
        T = (1+s1/ki)/(1+s0/ki); rec[f'T_Ki{ki}'] = round(T,2); rec[f'uptake_ratio_Ki{ki}'] = round(ra/T,2)
    u = rec['uptake_ratio_Ki300']
    rec['call'] = 'above threshold (uptake up)' if u > 1.1 else ('below threshold (uptake down)' if u < 0.9 else 'near threshold')
    g = r['obs_grain_pct']; st = r['obs_straw_pct']
    rec['obs_grain_pct'] = g; rec['obs_straw_pct'] = st
    sig = lambda x: '' if x == '' else ('+' if float(x) > 0 else '-')
    rec['grain_sign_matches'] = '' if g == '' else ('yes' if (u > 1) == (float(g) > 0) else 'no')
    rec['straw_sign_matches'] = '' if st == '' else ('yes' if (u > 1) == (float(st) > 0) else 'no')
    out.append(rec)
with open('threshold_table.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for o in out:
    print(f"{o['study'][:12]:12s} {o['treatment'][:26]:26s} Si x{o['si_ratio']:<5} As x{o['as_ratio']:<5} T {o['T_Ki300']:<5} U {o['uptake_ratio_Ki300']:<5} grain {o['obs_grain_pct']:>4} straw {o['obs_straw_pct']:>4}  {o['grain_sign_matches']}/{o['straw_sign_matches']}")
n = [o for o in out if o['grain_sign_matches']]; print('grain sign match', sum(o['grain_sign_matches']=='yes' for o in n), '/', len(n))
n = [o for o in out if o['straw_sign_matches']]; print('straw sign match', sum(o['straw_sign_matches']=='yes' for o in n), '/', len(n))
