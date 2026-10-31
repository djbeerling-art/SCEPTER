# Per-cell crop calendar and N fertiliser rate for the global riceland runs.
#   crop calendar  GGCMI phase 3 crop calendar (Jaegermeyr et al. 2021, Nature Food; zenodo 5062513), 0.5 degree,
#                  main rice season, irrigated (ri1_ir_*.nc4) and rainfed (ri1_rf_*.nc4): planting_day and
#                  growing_season_length (days). Built from RiceAtlas, MIRCA2000 and Sacks et al. 2010.
#   N fertiliser   EarthStat crop-specific fertiliser (Mueller et al. 2012, Nature), rice N application rate
#                  (kg N per ha of harvested rice, ~2000), 5 arcmin GeoTIFF (rice_NitrogenApplication_Rate.tif)
# Each cell takes the value of its own pixel, or of the nearest valid pixel within --fillkm (coastal and delta
# cells fall outside coarse land masks); cells still without an N rate take their country's area-weighted median,
# then --nfert0 (default 120 kg N/ha, the GMD runs' value).
# Usage: python3 paddy/global/cal_nfert.py <cells.csv> <ri1_ir.nc4> <ri1_rf.nc4> <rice_N_rate.tif> <out.csv>
#        [--fillkm=150] [--nfert0=120]
# Output: cell, plant_irr, gsl_irr, plant_rf, gsl_rf (day of year, days), nfert (kg N/ha per season), fill flags.
import sys
import numpy as np, pandas as pd

R_EARTH = 6371.

def opt(name, default):
    return next((a.split('=', 1)[1] for a in sys.argv if a.startswith(f'--{name}=')), default)

def nearest(lon_v, lat_v, val, lon_q, lat_q, fillkm):
    """Value at the nearest valid point (haversine) within fillkm, else NaN; also the distance in km."""
    from sklearn.neighbors import BallTree
    ok = np.isfinite(val)
    tree = BallTree(np.radians(np.c_[lat_v[ok], lon_v[ok]]), metric='haversine')
    d, i = tree.query(np.radians(np.c_[lat_q, lon_q]), k=1)
    d = d[:, 0] * R_EARTH
    out = val[ok][i[:, 0]].astype(float)
    out[d > fillkm] = np.nan
    return out, d

def grid_points(lon, lat, arr):
    LON, LAT = np.meshgrid(lon, lat)
    return LON.ravel(), LAT.ravel(), arr.ravel()

def read_calendar(path):
    import xarray as xr
    ds = xr.open_dataset(path, decode_times=False)
    lat = ds['lat'].values; lon = ds['lon'].values
    pd_ = ds['planting_day'].values.astype(float)
    gsl = ds['growing_season_length'].values.astype(float)
    pd_[(pd_ < 1) | (pd_ > 366)] = np.nan
    gsl[(gsl < 30) | (gsl > 300)] = np.nan
    return lon, lat, pd_, gsl

def read_tif(path):
    import rasterio
    with rasterio.open(path) as r:
        a = r.read(1).astype(float)
        if r.nodata is not None: a[a == r.nodata] = np.nan
        t = r.transform
        lon = t.c + t.a * (np.arange(r.width) + 0.5)
        lat = t.f + t.e * (np.arange(r.height) + 0.5)
    a[(a < 0) | (a > 1000)] = np.nan
    return lon, lat, a

def main(cells_f, ir_f, rf_f, n_f, out_f):
    fillkm = float(opt('fillkm', 150)); nfert0 = float(opt('nfert0', 120))
    c = pd.read_csv(cells_f)
    out = pd.DataFrame({'cell': c.cell})
    for tag, f in (('irr', ir_f), ('rf', rf_f)):
        lon, lat, p, g = read_calendar(f)
        for var, arr in (('plant', p), ('gsl', g)):
            lo, la, v = grid_points(lon, lat, arr)
            val, d = nearest(lo, la, v, c.lon.values, c.lat.values, fillkm)
            out[f'{var}_{tag}'] = np.round(val)
        out[f'cal_{tag}_km'] = np.round(d, 1)
    lo, la, v = grid_points(*read_tif(n_f))
    keep = np.isfinite(v) & (v > 0)          # zero = no rice in the source map, not an observed zero rate
    n, d = nearest(lo[keep], la[keep], v[keep], c.lon.values, c.lat.values, fillkm)
    out['nfert_src'] = np.where(np.isfinite(n), 'earthstat', '')
    area = c.area_irr + c.area_rf
    for iso, idx in c.groupby('iso3').groups.items():
        m = np.isfinite(n[idx])
        if m.any():
            vals = n[idx][m]; w = area.values[idx][m]
            o = np.argsort(vals); cw = np.cumsum(w[o]) / w.sum()
            med = vals[o][np.searchsorted(cw, 0.5)]
            fill = ~m
            n[np.asarray(idx)[fill]] = med
            out.loc[np.asarray(idx)[fill], 'nfert_src'] = 'country_median'
    miss = ~np.isfinite(n)
    n[miss] = nfert0; out.loc[miss, 'nfert_src'] = 'default'
    out['nfert'] = np.round(n, 1)
    # cells with no calendar at all (neither regime) fall back to the other regime, then to day 150 / 113 d
    for var, dflt in (('plant', 150), ('gsl', 113)):
        a, b = f'{var}_irr', f'{var}_rf'
        out[a] = out[a].fillna(out[b]); out[b] = out[b].fillna(out[a])
        out[[a, b]] = out[[a, b]].fillna(dflt)
    out.to_csv(out_f, index=False)
    w = area.values
    print(f'{len(out)} cells -> {out_f}')
    print('N source (area share):', {k: round(float(w[out.nfert_src == k].sum() / w.sum()), 3) for k in out.nfert_src.unique()})
    print(f'area-weighted N {np.average(out.nfert, weights=w):.0f} kg N/ha; '
          f'planting day irr median {out.plant_irr.median():.0f}, rf {out.plant_rf.median():.0f}')

if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    main(*a)
