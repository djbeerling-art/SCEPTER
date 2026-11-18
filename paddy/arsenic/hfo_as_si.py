"""PHREEQC prototype: As(III), As(V), Si and P competing for ferrihydrite (HFO) sites in a flooded paddy soil.

Dzombak & Morel (1990) diffuse-layer model with the official USGS wateq4f.dat (copied here from the phreeqc3
repository, May 2024 build). Its HFO constants include arsenate (29.31 / 23.51 / 10.58), arsenite (5.41),
phosphate (31.29 / 25.39 / 17.72) and silicic acid (4.28 / -3.22 / -11.69).
Dissolved Si and pH are imposed (fixed-activity phases), as SCEPTER would supply them. The bare surface is first
conditioned with a background electrolyte, so As and P totals are conserved between porewater and surfaces.

Base case per litre of porewater (all adjustable below):
  2.6 kg soil (bulk density 1.3, porosity 0.5); reactive Fe(III) 0.5 wt% as ferrihydrite -> 0.233 mol Fe/L
  sites: 0.2 mol weak + 0.005 mol strong per mol Fe, 600 m2/g, 89 g/mol
  labile As 6 mg/kg soil (of ~20 mg/kg total) -> 0.208 mmol/L, 95% as As(III) (fixed by low pe)
  labile P 1.0 mmol/L
"""
import sys, itertools, csv, math
from pathlib import Path
from phreeqpython import PhreeqPython
pp = PhreeqPython(database='wateq4f.dat', database_directory=Path(__file__).resolve().parent)
ip = pp.ip

EXTRA = """
SURFACE_MASTER_SPECIES
  Hfx_w Hfx_wOH
SURFACE_SPECIES
  # Si- and P-insensitive As pool (crystalline Fe, clays, OM); same As and acid-base constants as Hfo_w
  Hfx_wOH = Hfx_wOH ; log_k 0
  Hfx_wOH + H+ = Hfx_wOH2+ ; log_k 7.29
  Hfx_wOH = Hfx_wO- + H+ ; log_k -8.93
  Hfx_wOH + H3AsO3 = Hfx_wH2AsO3 + H2O ; log_k 5.41
  Hfx_wOH + AsO4-3 + 3H+ = Hfx_wH2AsO4 + H2O ; log_k 29.31
  Hfx_wOH + AsO4-3 + 2H+ = Hfx_wHAsO4- + H2O ; log_k 23.51
  Hfx_wOH + AsO4-3 = Hfx_wOHAsO4-3 ; log_k 10.58
PHASES
  Fix_H+
  H+ = H+ ; log_k 0
  Fix_Si
  H4SiO4 = H4SiO4 ; log_k 0
"""

def run(pH=6.8, si_mM=0.3, p_mM=1.0, as_mM=0.208, fe_mol=0.233, as5_frac=0.05, xs=0.0):
    """xs = Si-insensitive sites as a multiple of ferrihydrite weak sites (0 = ferrihydrite only)."""
    w, s = 0.2 * fe_mol, 0.005 * fe_mol
    mass = fe_mol * 89.0
    inp = f"""{EXTRA}
SOLUTION 1  # background electrolyte used to condition the bare surface
  units mmol/kgw
  pH 6.5
  Na 5
  Cl 5 charge
  Ca 2
SURFACE 1
  Hfo_sOH {s:.6g} 600 {mass:.6g}
  Hfo_wOH {w:.6g}
{f'  Hfx_wOH {xs*w:.6g} 600 {xs*mass:.6g}' if xs > 0 else ''}
  -equilibrate 1
SOLUTION 2  # porewater carrying the labile As and P; As and P totals are conserved between solution and surfaces
  units mmol/kgw
  pH 6.5
  Na 5
  Cl 5 charge
  Ca 2
  As(3) {as_mM*(1-as5_frac):.6g}
  As(5) {as_mM*as5_frac:.6g}
  P {p_mM:.6g}
  Si {si_mM:.6g}
EQUILIBRIUM_PHASES 1
  Fix_H+ {-pH} {'HCl' if pH < 6.5 else 'NaOH'} 10
  Fix_Si {math.log10(si_mM*1e-3):.4f} SiO2 10
END
USE solution 2
USE surface 1
USE equilibrium_phases 1
SELECTED_OUTPUT
  -reset false
  -pH true
USER_PUNCH
  -headings As3_uM As5_uM P_uM Si_uM
  10 PUNCH TOT("As(3)")*1e6, TOT("As(5)")*1e6, TOT("P")*1e6, TOT("Si")*1e6
END
"""
    try:
        ip.run_string(inp)
    except Exception:
        return None
    out = ip.get_selected_output_array()
    d = dict(zip(out[0], out[-1]))
    if abs(d['pH'] - pH) > 0.02 or abs(d['Si_uM'] / (si_mM * 1e3) - 1) > 0.15:
        return None
    d['As_sorbed_frac'] = 1 - (d['As3_uM'] + d['As5_uM']) / (as_mM * 1e3)
    return d

def fmt(r, k, n):
    return '' if r is None else round(r[k], n)

if __name__ == '__main__':
    w = csv.writer(sys.stdout)
    w.writerow(['scan', 'pH', 'Si_mM', 'P_mM', 'Fe_frac', 'As3_uM', 'As5_uM', 'P_uM', 'As_sorbed_frac'])
    for si in (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5):
        for pH in (6.0, 6.8, 7.5):
            r = run(pH=pH, si_mM=si); w.writerow(['Si', pH, si, 1.0, 1.0, fmt(r, 'As3_uM', 3), fmt(r, 'As5_uM', 4), fmt(r, 'P_uM', 2), fmt(r, 'As_sorbed_frac', 4)])
    for pH in (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0):
        r = run(pH=pH); w.writerow(['pH', pH, 0.3, 1.0, 1.0, fmt(r, 'As3_uM', 3), fmt(r, 'As5_uM', 4), fmt(r, 'P_uM', 2), fmt(r, 'As_sorbed_frac', 4)])
    for p in (0.1, 0.3, 1.0, 3.0, 10.0):
        r = run(p_mM=p); w.writerow(['P', 6.8, 0.3, p, 1.0, fmt(r, 'As3_uM', 3), fmt(r, 'As5_uM', 4), fmt(r, 'P_uM', 2), fmt(r, 'As_sorbed_frac', 4)])
    for f in (1.0, 0.8, 0.6, 0.4, 0.2, 0.1):
        r = run(fe_mol=0.233 * f); w.writerow(['Fe_loss', 6.8, 0.3, 1.0, f, fmt(r, 'As3_uM', 3), fmt(r, 'As5_uM', 4), fmt(r, 'P_uM', 2), fmt(r, 'As_sorbed_frac', 4)])
