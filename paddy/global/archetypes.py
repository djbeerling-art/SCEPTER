# Group rice cells into model archetypes and map each archetype to SCEPTER driver settings.
# Input: cells CSV, one row per rice cell (e.g. 5 arcmin), columns
#   cell, lon, lat, iso3, area_irr, area_rf (ha, MapSPAM), ph (H2O, 0-15 cm), cec (cmol/kg), soc (wt %),
#   temp (growing-season mean, deg C), rock (name of a data/dust_<rock>.in recipe, from basalt_norm.py for the
#   nearest basalt outcrop within 250 km); optional clay (wt %, SoilGrids) or fe3 (reducible Fe(III) wt fraction)
# Soil reducible Fe(III) (third-ranked control on dCH4 in the Morris screen) has no global map: if fe3 is not given it
# is scaled from clay, fe3 = 0.012*clay/30 clipped to 0.004-0.03 (a placeholder proxy; run global results at the
# clip bounds too for its uncertainty), else the reference 0.012.
# Each cell splits into an irrigated part (continuous flooding with a mid-season drain, --phd=14) and a rainfed
# part (--awd as a stand-in for intermittent flooding). Within each (regime, rock) group, cells are clustered by
# k-means on standardised (ph, cec, soc, temp[, fe3]), weighted by rice area; archetypes are allotted to groups in
# proportion to area (at least one each). Archetype values are area-weighted means of member cells.
# Usage: python3 archetypes.py cells.csv outdir [--n=200]   (or <archetypes.csv> --reargs <out.csv>, see below)
#   -> outdir/archetypes.csv (one row per archetype, with driver settings) and outdir/cell_map.csv
import sys, os
import numpy as np, pandas as pd
from sklearn.cluster import KMeans

FEAT = ['ph', 'cec', 'soc', 'temp']
# Soil-conditioning rule of the v22 GMD archetype runs (29 Oct 2026): dry-soil pH >= 6.8 -> calcareous soil
# (--soilcc=0.005 --irralk=10); otherwise irrigation alkalinity interpolated on pH [5.0, 6.0, 6.2, 6.5] -> [0, 1, 4, 10].
# (Replaces the v5 15-cell spin-up lookup.) Reproduces the GMD args, e.g. irr061 pH 6.04 -> 1.6, rf189 pH 6.34 -> 6.8.
PH_IRRALK = [(5.0, 0.), (6.0, 1.), (6.2, 4.), (6.5, 10.)]   # (pH, irralk)
PH_CALC, SOILCC_CALC, IRRALK_CALC = 6.8, 0.005, 10.
SOC_REF, OM_REF = 1.5, 450.                         # wt % SOC of the reference soil; its OM input (g C/m2/yr)
FE3_REF, CLAY_REF = 0.012, 30.                      # reference soil reducible Fe(III) (wt fraction) and assumed clay %

def ph_to_soil(ph):
    """Dry-soil pH -> (irralk, soilcc), v22 rule; clamped at both ends."""
    if ph >= PH_CALC: return IRRALK_CALC, SOILCC_CALC
    y, x = zip(*PH_IRRALK)
    return float(np.interp(ph, y, x)), 0.

def driver_args(r):
    irralk, soilcc = ph_to_soil(r.ph)
    om = OM_REF*float(np.clip(r.soc/SOC_REF, 0.5, 2.0))
    a = f'--rock={r.rock} --irralk={irralk:.2f} --cec={r.cec:.1f} --temp={r.temp:.1f} --om={om:.0f} --fe3={r.fe3:.4f}'
    if soilcc > 0: a += f' --soilcc={soilcc:.3f}'
    return a

def build(cells, n):
    cells = cells.copy()
    if 'fe3' not in cells:
        cells['fe3'] = np.clip(FE3_REF*cells.clay/CLAY_REF, 0.004, 0.03) if 'clay' in cells else FE3_REF
    feat = FEAT+['fe3'] if cells.fe3.nunique() > 1 else FEAT
    rows = []
    for regime, col in (('irr', 'area_irr'), ('rf', 'area_rf')):
        c = cells[cells[col] > 0].copy(); c['area'] = c[col]; c['regime'] = regime
        rows.append(c)
    c = pd.concat(rows, ignore_index=True)
    mu, sd = c[feat].mean(), c[feat].std().replace(0, 1)
    groups = c.groupby(['regime', 'rock'])
    tot = c.area.sum()
    arch, cmap = [], []
    for (regime, rock), g in groups:
        k = int(min(len(g), max(1, round(n*g.area.sum()/tot))))
        lab = KMeans(k, n_init=10, random_state=0).fit(((g[feat]-mu)/sd).values, sample_weight=g.area.values).labels_
        for j in range(k):
            m = g[lab == j]; w = m.area.values
            aid = f'{regime}{len(arch):03d}'
            arch.append(dict(arch=aid, regime=regime, rock=rock, area=w.sum(), ncell=len(m),
                             **{f: np.average(m[f], weights=w) for f in FEAT+['fe3']}))
            cmap += [(cl, regime, aid, a) for cl, a in zip(m.cell, w)]
    a = pd.DataFrame(arch)
    a['args'] = [driver_args(r) for r in a.itertuples()]
    return a, pd.DataFrame(cmap, columns=['cell', 'regime', 'arch', 'area'])

if __name__ == '__main__' and '--reargs' in sys.argv:
    # Recompute only the driver args of an existing archetypes.csv (keeps archetype ids/membership):
    #   python3 archetypes.py <archetypes.csv in> --reargs <archetypes.csv out>
    a = pd.read_csv(sys.argv[1]); a['args'] = [driver_args(r) for r in a.itertuples()]
    a.to_csv(sys.argv[3], index=False); print(f'{len(a)} archetypes: args recomputed'); sys.exit()

if __name__ == '__main__':
    cells = pd.read_csv(sys.argv[1]); out = sys.argv[2]
    n = int(next((x.split('=')[1] for x in sys.argv if x.startswith('--n=')), 200))
    os.makedirs(out, exist_ok=True)
    a, m = build(cells, n)
    a.to_csv(os.path.join(out, 'archetypes.csv'), index=False); m.to_csv(os.path.join(out, 'cell_map.csv'), index=False)
    print(f'{len(cells)} cells -> {len(a)} archetypes ({(a.regime=="irr").sum()} irrigated, {(a.regime=="rf").sum()} rainfed)')
