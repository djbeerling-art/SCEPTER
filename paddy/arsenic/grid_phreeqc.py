"""Generate a PHREEQC (Dzombak & Morel DLM) grid of As(III), As(V) and Si sorption on ferrihydrite
for fitting SCEPTER's non-electrostatic competitive scheme. P held low (0.01 mM; SCEPTER carries no P).
Usage: python3 grid_phreeqc.py [P] -> grid_phreeqc.csv (P 0.01 mM) or grid_phreeqc_P<P>.csv.
P is total labile P in mM per L porewater, or 'r<x>' for P = x mol per mol ferrihydrite Fe (labile P scales with
reactive Fe in soils; r0.02 and r0.1 bracket low and high P loading)."""
import sys, csv, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hfo_as_si import run
PA = sys.argv[1] if len(sys.argv) > 1 else '0.01'
ratio = PA.startswith('r'); P = float(PA[1:] if ratio else PA)
rows = []
for pH, si, asm, fe, f5 in itertools.product(
        (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0),
        (0.02, 0.1, 0.3, 0.6, 1.0, 1.5, 2.0),
        (0.002, 0.02, 0.2),
        (0.233, 0.07, 0.023),
        (0.05, 0.5)):
    r = run(pH=pH, si_mM=si, p_mM=(P*fe*1e3 if ratio else P), as_mM=asm, fe_mol=fe, as5_frac=f5)
    if r is None: continue
    rows.append(dict(pH=pH, si_mM=si, as_mM=asm, fe_mol=fe, as5_frac=f5,
                     as3_uM=r['As3_uM'], as5_uM=r['As5_uM'], si_uM=r['Si_uM']))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'grid_phreeqc.csv' if PA == '0.01' else f'grid_phreeqc_P{PA}.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows))
