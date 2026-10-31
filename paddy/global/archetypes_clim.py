# Climate-aware archetypes for the monthly-forcing global runs (successor of archetypes.py, which stays as the
# frozen v22/GMD grouping).
# Inputs
#   cells.csv      data/global/cells.csv (soil, rock, rice area; see archetypes.py)
#   cal.csv        paddy/global/cal_nfert.py output: plant_irr, plant_rf (day of year), nfert (kg N/ha per season)
#   forcing        monthly per-cell forcing, CSV(.gz) or parquet, long format: cell, month, tsoil_C, wsat, p_mm, et_mm
#                  (era5land_clim_cells.csv, or forcing_cells.parquet with gcm/ssp/period columns: the baseline rows
#                  gcm=ERA5-Land are used for clustering)
# Clustering: within each (regime, rock) group, area-weighted k-means on standardised
#   soil   ph, cec, soc, fe3 (clay proxy as in archetypes.py), nfert (capped at --nmax, default 300)
#   timing sin and cos of the planting day (so day 360 and day 5 are neighbours)
#   climate mean soil T over the 113-d season from planting, mean soil T, wetness (wsat) and P-ET over the fallow
# Archetype values are area-weighted means of member cells (planting day as a circular mean). Each archetype gets
#   <outdir>/clim/<case>/<arch>.csv  (month, tsoil_C, wsat, p_mm, et_mm: area-weighted monthly means of its cells)
#   and driver args with --plant, --nfert and --clim=<that file> added (--temp = season-mean soil T, for reference).
# The 113-d season length is kept for every archetype (David, 31 Oct 2026); GGCMI sets only the planting day.
# Usage
#   python3 archetypes_clim.py <cells.csv> <cal.csv> <forcing> <outdir> [--n=600] [--nmax=300]
#       -> outdir/archetypes.csv, outdir/cell_map.csv, outdir/clim/baseline/*.csv
#   python3 archetypes_clim.py --case <outdir> <forcing_cells.parquet> <gcm> <ssp> <period>
#       -> outdir/clim/<gcm>_<ssp>_<period>/*.csv and outdir/archetypes_<gcm>_<ssp>_<period>.csv (same archetypes
#          and membership; --clim points at the case's files and --temp is that case's season mean)
import sys, os
import numpy as np, pandas as pd
from sklearn.cluster import KMeans
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archetypes import ph_to_soil, OM_REF, SOC_REF, FE3_REF, CLAY_REF

SEASON_D = 113
MID = 15.2 + 30.42*np.arange(12)            # mid-month day of year
CLIM_COLS = ['tsoil_C', 'wsat', 'p_mm', 'et_mm']

def opt(name, default):
    return next((a.split('=', 1)[1] for a in sys.argv if a.startswith(f'--{name}=')), default)

def read_forcing(path, case=None):
    cols = ['cell', 'month'] + CLIM_COLS
    if path.endswith('.parquet'):
        import pyarrow.parquet as pq
        names = pq.read_schema(path).names
        filt = None
        if 'gcm' in names:
            g, s, p = case or ('ERA5-Land', 'historical', '1995-2014')
            filt = [('gcm', '=', g), ('ssp', '=', s), ('period', '=', p)]
        f = pd.read_parquet(path, columns=cols, filters=filt)
    else:
        f = pd.read_csv(path, usecols=cols)
    f = f.sort_values(['cell', 'month'])
    assert (f.groupby('cell').size() == 12).all(), 'forcing needs 12 months per cell'
    return f

def monthly(f):
    """cell index and (ncell, 12) arrays per forcing column."""
    cells = f.cell.values.reshape(-1, 12)[:, 0]
    return cells, {c: f[c].values.reshape(-1, 12).astype(float) for c in CLIM_COLS}

def window_mean(arr, start, length):
    """Mean of each row's periodic monthly series (linear between mid-months) over [start, start+length) days."""
    days = (start[:, None] + np.arange(0, length, 5)[None, :]) % 365
    ext_d = np.r_[MID - 365, MID, MID + 365]
    out = np.empty(len(arr))
    for i in range(len(arr)):
        out[i] = np.interp(days[i], ext_d, np.tile(arr[i], 3)).mean()
    return out

def climate_features(plant, M):
    fallow0 = (plant + SEASON_D) % 365
    q = np.maximum(M['p_mm'] - M['et_mm'], 0)
    return pd.DataFrame({
        't_season': window_mean(M['tsoil_C'], plant, SEASON_D),
        't_fallow': window_mean(M['tsoil_C'], fallow0, 365 - SEASON_D),
        'w_fallow': window_mean(M['wsat'], fallow0, 365 - SEASON_D),
        'q_fallow': window_mean(q, fallow0, 365 - SEASON_D)})

def circ_mean_day(days, w):
    a = 2*np.pi*np.asarray(days)/365
    return (np.degrees(np.arctan2(np.average(np.sin(a), weights=w), np.average(np.cos(a), weights=w))) % 360)/360*365

def driver_args(r, clim_path):
    irralk, soilcc = ph_to_soil(r.ph)
    om = OM_REF*float(np.clip(r.soc/SOC_REF, 0.5, 2.0))
    a = (f'--rock={r.rock} --irralk={irralk:.2f} --cec={r.cec:.1f} --temp={r.t_season:.1f} --om={om:.0f} '
         f'--fe3={r.fe3:.4f} --nfert={r.nfert:.0f} --plant={r.plant:.0f} --clim={clim_path}')
    if soilcc > 0: a += f' --soilcc={soilcc:.3f}'
    return a

def write_clim(cmap, cells_idx, M, outdir, case):
    """Area-weighted monthly forcing per archetype; returns {arch: path} and per-archetype season-mean T inputs."""
    d = os.path.join(outdir, 'clim', case); os.makedirs(d, exist_ok=True)
    pos = pd.Series(np.arange(len(cells_idx)), index=cells_idx)
    paths, mon = {}, {}
    for aid, g in cmap.groupby('arch'):
        i = pos.loc[g.cell].values; w = g.area.values
        df = pd.DataFrame({'month': np.arange(1, 13), **{c: np.average(M[c][i], axis=0, weights=w) for c in CLIM_COLS}})
        p = os.path.abspath(os.path.join(d, f'{aid}.csv'))
        df.round(4).to_csv(p, index=False); paths[aid] = p; mon[aid] = df
    return paths, mon

def build(cells, cal, forc, outdir, n, nmax):
    c = cells.merge(cal[['cell', 'plant_irr', 'plant_rf', 'nfert']], on='cell', how='left')
    c['nfert'] = c.nfert.fillna(120).clip(upper=nmax)
    if 'fe3' not in c:
        c['fe3'] = np.clip(FE3_REF*c.clay/CLAY_REF, 0.004, 0.03) if 'clay' in c else FE3_REF
    fcells, M = monthly(forc)
    bad = np.zeros(len(fcells), bool)
    for v in M.values(): bad |= ~np.isfinite(v).all(axis=1)
    have = set(fcells[~bad])
    miss = ~c.cell.isin(have)
    if miss.any():
        print(f'WARNING: {miss.sum()} cells ({c.loc[miss, ["area_irr", "area_rf"]].sum().sum()/1e6:.2f} Mha) have no (complete) forcing; dropped')
        c = c[~miss]
    pos = pd.Series(np.arange(len(fcells)), index=fcells)
    rows = []
    for regime, col in (('irr', 'area_irr'), ('rf', 'area_rf')):
        g = c[c[col] > 0].copy(); g['area'] = g[col]; g['regime'] = regime
        g['plant'] = g[f'plant_{regime}'].fillna(150).values
        idx = pos.loc[g.cell].values
        g = pd.concat([g.reset_index(drop=True),
                       climate_features(g.plant.values, {k: v[idx] for k, v in M.items()})], axis=1)
        rows.append(g)
    c = pd.concat(rows, ignore_index=True)
    c['psin'] = np.sin(2*np.pi*c.plant/365); c['pcos'] = np.cos(2*np.pi*c.plant/365)
    feat = ['ph', 'cec', 'soc', 'fe3', 'nfert', 'psin', 'pcos', 't_season', 't_fallow', 'w_fallow', 'q_fallow']
    mu, sd = c[feat].mean(), c[feat].std().replace(0, 1)
    tot = c.area.sum(); arch, cmap = [], []
    for (regime, rock), g in c.groupby(['regime', 'rock']):
        k = int(min(len(g), max(1, round(n*g.area.sum()/tot))))
        lab = KMeans(k, n_init=4, random_state=0).fit(((g[feat]-mu)/sd).values, sample_weight=g.area.values).labels_
        for j in range(k):
            m = g[lab == j]; w = m.area.values
            aid = f'{regime}{len(arch):04d}'
            rec = dict(arch=aid, regime=regime, rock=rock, area=w.sum(), ncell=len(m),
                       plant=circ_mean_day(m.plant, w),
                       **{f: np.average(m[f], weights=w) for f in ['ph', 'cec', 'soc', 'fe3', 'nfert', 'clay',
                                                                   't_season', 't_fallow', 'w_fallow', 'q_fallow']
                          if f in m})
            arch.append(rec); cmap += [(cl, regime, aid, a) for cl, a in zip(m.cell, w)]
    a = pd.DataFrame(arch); cm = pd.DataFrame(cmap, columns=['cell', 'regime', 'arch', 'area'])
    paths, _ = write_clim(cm, fcells, M, outdir, 'baseline')
    a['args'] = [driver_args(r, paths[r.arch]) for r in a.itertuples()]
    return a, cm

def case_main(outdir, forcing, gcm, ssp, period):
    a = pd.read_csv(os.path.join(outdir, 'archetypes.csv')); cm = pd.read_csv(os.path.join(outdir, 'cell_map.csv'))
    fcells, M = monthly(read_forcing(forcing, (gcm, ssp, period)))
    case = f'{gcm}_{ssp}_{period}'
    paths, mon = write_clim(cm, fcells, M, outdir, case)
    # season-mean soil T of the case for --temp (reference only; --clim drives the run)
    a['t_season'] = [window_mean(mon[r.arch].tsoil_C.values[None, :], np.array([r.plant]), SEASON_D)[0]
                     for r in a.itertuples()]
    a['args'] = [driver_args(r, paths[r.arch]) for r in a.itertuples()]
    out = os.path.join(outdir, f'archetypes_{case}.csv'); a.to_csv(out, index=False)
    print(f'{case}: {len(a)} archetypes -> {out}')

if __name__ == '__main__':
    if '--case' in sys.argv:
        x = [v for v in sys.argv[1:] if not v.startswith('--')]
        case_main(*x); sys.exit()
    x = [v for v in sys.argv[1:] if not v.startswith('--')]
    cells_f, cal_f, forc_f, out = x
    n = int(opt('n', 600)); nmax = float(opt('nmax', 300))
    os.makedirs(out, exist_ok=True)
    a, m = build(pd.read_csv(cells_f), pd.read_csv(cal_f), read_forcing(forc_f), out, n, nmax)
    a.to_csv(os.path.join(out, 'archetypes.csv'), index=False); m.to_csv(os.path.join(out, 'cell_map.csv'), index=False)
    w = a.area
    print(f'{m.cell.nunique()} cells -> {len(a)} archetypes ({(a.regime == "irr").sum()} irrigated, '
          f'{(a.regime == "rf").sum()} rainfed); area-weighted season soil T {np.average(a.t_season, weights=w):.1f} C, '
          f'N {np.average(a.nfert, weights=w):.0f} kg/ha')
