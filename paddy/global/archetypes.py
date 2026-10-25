# Group rice cells into model archetypes and map each archetype to SCEPTER driver settings.
# Input: cells CSV, one row per rice cell (e.g. 5 arcmin), columns
#   cell, lon, lat, iso3, area_irr, area_rf (ha, MapSPAM), ph (H2O, 0-15 cm), cec (cmol/kg), soc (wt %),
#   temp (growing-season mean, deg C), rock (name of a data/dust_<rock>.in recipe, from basalt_norm.py for the
#   nearest basalt outcrop within 250 km)
# Each cell splits into an irrigated part (continuous flooding with a mid-season drain, --phd=14) and a rainfed
# part (--awd as a stand-in for intermittent flooding). Within each (regime, rock) group, cells are clustered by
# k-means on standardised (ph, cec, soc, temp), weighted by rice area; archetypes are allotted to groups in
# proportion to area (at least one each). Archetype values are area-weighted means of member cells.
# Usage: python3 archetypes.py cells.csv outdir [--n=200]
#   -> outdir/archetypes.csv (one row per archetype, with driver settings) and outdir/cell_map.csv
import sys, os
import numpy as np, pandas as pd
from sklearn.cluster import KMeans

FEAT = ['ph', 'cec', 'soc', 'temp']
# spin-up topsoil pH (3-yr fallow, 0-15 cm median) of the reference soil against irrigation alkalinity and soil calcite
PH_IRRALK = [(0.0, 4.95), (0.5, 6.38), (1.0, 6.57), (2.0, 6.72)]
PH_SOILCC = [(0.0, 6.72), (0.027, 7.45)]           # with irralk 2
SOC_REF, OM_REF = 1.5, 450.                         # wt % SOC of the reference soil; its OM input (g C/m2/yr)

def ph_to_soil(ph):
    """Target pH -> (irralk, soilcc) by linear interpolation of the spin-up lookup; clamped at both ends."""
    if ph <= PH_IRRALK[-1][1]:
        x, y = zip(*PH_IRRALK)
        return float(np.interp(ph, y, x)), 0.
    x, y = zip(*PH_SOILCC)
    return PH_IRRALK[-1][0], float(np.interp(ph, y, x))

def driver_args(r):
    irralk, soilcc = ph_to_soil(r.ph)
    om = OM_REF*float(np.clip(r.soc/SOC_REF, 0.5, 2.0))
    a = f'--rock={r.rock} --irralk={irralk:.3f} --cec={r.cec:.1f} --temp={r.temp:.1f} --om={om:.0f}'
    if soilcc > 0: a += f' --soilcc={soilcc:.4f}'
    return a

def build(cells, n):
    rows = []
    for regime, col in (('irr', 'area_irr'), ('rf', 'area_rf')):
        c = cells[cells[col] > 0].copy(); c['area'] = c[col]; c['regime'] = regime
        rows.append(c)
    c = pd.concat(rows, ignore_index=True)
    mu, sd = c[FEAT].mean(), c[FEAT].std().replace(0, 1)
    groups = c.groupby(['regime', 'rock'])
    tot = c.area.sum()
    arch, cmap = [], []
    for (regime, rock), g in groups:
        k = int(min(len(g), max(1, round(n*g.area.sum()/tot))))
        lab = KMeans(k, n_init=10, random_state=0).fit(((g[FEAT]-mu)/sd).values, sample_weight=g.area.values).labels_
        for j in range(k):
            m = g[lab == j]; w = m.area.values
            aid = f'{regime}{len(arch):03d}'
            arch.append(dict(arch=aid, regime=regime, rock=rock, area=w.sum(), ncell=len(m),
                             **{f: np.average(m[f], weights=w) for f in FEAT}))
            cmap += [(cl, regime, aid, a) for cl, a in zip(m.cell, w)]
    a = pd.DataFrame(arch)
    a['args'] = [driver_args(r) for r in a.itertuples()]
    return a, pd.DataFrame(cmap, columns=['cell', 'regime', 'arch', 'area'])

if __name__ == '__main__':
    cells = pd.read_csv(sys.argv[1]); out = sys.argv[2]
    n = int(next((x.split('=')[1] for x in sys.argv if x.startswith('--n=')), 200))
    os.makedirs(out, exist_ok=True)
    a, m = build(cells, n)
    a.to_csv(os.path.join(out, 'archetypes.csv'), index=False); m.to_csv(os.path.join(out, 'cell_map.csv'), index=False)
    print(f'{len(cells)} cells -> {len(a)} archetypes ({(a.regime=="irr").sum()} irrigated, {(a.regime=="rf").sum()} rainfed)')
