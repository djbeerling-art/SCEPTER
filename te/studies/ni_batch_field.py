"""Does equilibrium chemistry alone explain the field decline in exchangeable Ni?

Closed batch per kg soil (0.3 kg water), soil-like solution at pCO2 = 10^-2,
Gaines-Thomas exchange (CEC from the plan), D&M Hfo sites from the reactive
fraction of soil Fe oxide. Ni total is fixed so that the control state has the
measured NH4OAc Ni (taken here as exchange + dissolved). We then raise pH
(and exchangeable Ca, as basalt does) and see where the Ni goes.
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from chem import thermo, system, solver, phreeqc_ref as pr

db = thermo.read([pr.DB, pr.EXTRA])
ELS = ['C', 'Ca', 'Cl', 'K', 'Mg', 'Na', 'Ni']
sy = system.build(db, ELS, valences={'C': ['C(4)']}, exchangers=['X'],
                  surfaces=['Hfo_w', 'Hfo_s'], phases=['CO2(g)'])
sv = solver.Speciator(sy)
W = 0.3
MNI = 58.69


def _surf(fe):
    return {'Hfo': solver.Surface({'Hfo_w': 0.2 * fe, 'Hfo_s': 0.005 * fe}, 600, 89 * fe)}


def state(pH, ni, cec, fe_mol, ca_mM, react):
    """Soil solution (charge-balanced, Ca-Mg-K-Na-Cl at pCO2 10^-2) sets the
    exchanger composition, as in a PHREEQC '-equilibrate'; Ni is then added
    as a fixed reactive total and partitioned."""
    fe = fe_mol * react
    ca = ca_mM * 1e-3
    sol = {'Ca': ca, 'Mg': ca / 4, 'K': 1e-4, 'Na': 2e-4}
    sol['Cl'] = 2 * ca + 2 * ca / 4 + 1e-4 + 2e-4
    # 1) exchanger in equilibrium with that solution (huge water mass)
    big = 1e6
    r0 = sv.solve({k: v * big for k, v in sol.items()}, pH=pH, W=big, cec={'X': cec},
                  fix_si=[('CO2(g)', -2.0, 'CO3-2')])
    tot = {}
    for k, v in sol.items():
        sp = sy.aliases[k]
        j = sy.index(sp)
        ex = sum(n * sy.ex.nu[sy.ex.names.index(nm), j] for nm, n in r0.ex.items())
        tot[k] = v * W + ex
    tot['Ni'] = ni
    r = sv.solve(tot, pH=pH, T=298.15, W=W, cec={'X': cec}, surfaces=_surf(fe),
                 fix_si=[('CO2(g)', -2.0, 'CO3-2')])
    return r.partition['Ni'], r


def calib(pH, target_mgkg, **kw):
    lo, hi = 1e-8, 1e-1
    for _ in range(60):
        mid = (lo * hi) ** 0.5
        p, _ = state(pH, mid, **kw)
        if (p['exchange'] + p['dissolved']) * MNI * 1e3 > target_mgkg:
            hi = mid
        else:
            lo = mid
    return mid


def run(label, cec, fe_oxide_pct, ph0, ni0, ph_list, ca0, ca1, reacts):
    fe_mol = fe_oxide_pct * 10 / 159.7 * 2   # mol Fe per kg soil (as Fe2O3)
    print(f'\n== {label}: CEC {cec*100:.0f} cmol/kg, Fe oxide {fe_oxide_pct}% '
          f'({fe_mol:.2f} mol Fe/kg), control pH {ph0}, NH4OAc Ni {ni0} mg/kg')
    for react in reacts:
        kw = dict(cec=cec, fe_mol=fe_mol, react=react)
        ni = calib(ph0, ni0, ca_mM=ca0, **kw)
        p0, r0 = state(ph0, ni, ca_mM=ca0, **kw)
        print(f'  reactive Fe fraction {react}: reactive Ni {ni*MNI*1e3:.1f} mg/kg, '
              f'control split exch {p0["exchange"]/ni:.2f} sorbed {p0["sorbed"]/ni:.2f} '
              f'diss {p0["dissolved"]/ni:.3f} (porewater {p0["dissolved"]/W*MNI*1e6:.0f} ug/L)')
        for ph in ph_list:
            for ca in (ca0, ca1):
                p, r = state(ph, ni, ca_mM=ca, **kw)
                lab = (p['exchange'] + p['dissolved']) * MNI * 1e3
                print(f'    pH {ph:.2f} Ca {ca:.1f} mM: NH4OAc-like Ni {lab:5.2f} mg/kg '
                      f'({lab/ni0:4.2f}x), sorbed {p["sorbed"]/ni:.2f}, '
                      f'porewater {p["dissolved"]/W*MNI*1e6:6.1f} ug/L')


if __name__ == '__main__':
    # Sabah P2 Acrisol: B plots 2018 pH 5.51, Ni 2.34 -> 2023 pH 6.01 Ni 0.89, 2024 pH 5.74 Ni 0.53
    run('Sabah (P2)', cec=0.06, fe_oxide_pct=3, ph0=5.5, ni0=2.34,
        ph_list=[5.5, 5.75, 6.0], ca0=0.5, ca1=2.0, reacts=[0.03, 0.1, 0.3])
    # Corn Belt T2 Mollisol: NH4OAc Ni ~0.09 mg/kg falling to ~0.05; pH ~6.2 rising ~0.2-0.4
    run('Corn Belt (T2)', cec=0.20, fe_oxide_pct=1.5, ph0=6.2, ni0=0.09,
        ph_list=[6.2, 6.4, 6.6], ca0=1.0, ca1=3.0, reacts=[0.03, 0.1, 0.3])
