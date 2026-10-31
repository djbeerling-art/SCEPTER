"""Monthly CMIP6 soil-temperature change (tsl) at rice cells.

For each GCM x SSP x period: delta(month) = mean tsl over the future period minus the same model's
historical 1995-2014 mean, at 0.10 m depth (linear interpolation between model soil layers), in K.
Ocean/non-land GCM cells are NaN; each rice cell takes the inverse-distance weighted mean of its
4 nearest GCM land points. Source: Pangeo CMIP6 Google Cloud archive (gs://cmip6, anonymous).

Usage: python3 cmip6_tsl_delta.py <pangeo-cmip6.csv> <cells.csv> <out_dir>
Writes tsl_delta_cells.nc (cell x gcm x ssp x period x month) and grid_<gcm>.nc (deltas on each GCM grid).
"""
import sys, os
import numpy as np, pandas as pd, xarray as xr
from scipy.spatial import cKDTree

GCMS = {'GFDL-ESM4': 'r1i1p1f1', 'IPSL-CM6A-LR': 'r1i1p1f1', 'MIROC6': 'r1i1p1f1',
        'BCC-CSM2-MR': 'r1i1p1f1', 'UKESM1-0-LL': 'r1i1p1f2'}
SSPS = ['ssp245', 'ssp585']
PERIODS = {'2041-2060': (2041, 2060), '2061-2080': (2061, 2080)}
HIST = (1995, 2014)
ZT = 0.10  # m

def xyz(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]

def open_tsl(cat, gcm, exp, mem):
    s = cat[(cat.variable_id == 'tsl') & (cat.table_id == 'Lmon') & (cat.source_id == gcm)
            & (cat.experiment_id == exp) & (cat.member_id == mem)]
    z = s.sort_values('version').zstore.iloc[-1]
    print(' ', z, flush=True)
    v = xr.open_zarr(z, storage_options={'token': 'anon'}, consolidated=True).tsl
    dz = [d for d in v.dims if d not in ('time', 'lat', 'lon')][0]
    return v.rename({dz: 'depth'})

def clim(v, y0, y1):
    yr = v.time.dt.year
    v = v.isel(time=((yr >= y0) & (yr <= y1)).values)   # works for 360-day calendars too
    assert v.time.size == 12 * (y1 - y0 + 1), (v.time.size, y0, y1)
    d = v.depth.values
    k = max(1, int(np.searchsorted(d, ZT)))          # layer below ZT (or 2nd layer if ZT above the top one)
    w = np.clip((ZT - d[k - 1]) / (d[k] - d[k - 1]), 0, 1)
    v = v.isel(depth=[k - 1, k]).load()
    v = v.isel(depth=0) * (1 - w) + v.isel(depth=1) * w
    return v.groupby('time.month').mean('time')

def main(catf, cellf, out):
    os.makedirs(out, exist_ok=True)
    cat = pd.read_csv(catf)
    cells = pd.read_csv(cellf, usecols=['cell', 'lon', 'lat'])
    q = xyz(cells.lon.values, cells.lat.values)
    res = np.full((len(cells), len(GCMS), len(SSPS), len(PERIODS), 12), np.nan, np.float32)
    for gi, (gcm, mem) in enumerate(GCMS.items()):
        print(gcm, flush=True)
        gf = f'{out}/grid_{gcm}.nc'
        if os.path.exists(gf):                      # resume: reuse deltas already computed for this GCM
            gds = xr.open_dataset(gf, engine='h5netcdf').load()
            grids = {k: gds[k] for k in gds.data_vars}
        else:
            h = clim(open_tsl(cat, gcm, 'historical', mem), *HIST)
            grids = {}
            for ssp in SSPS:
                f = open_tsl(cat, gcm, ssp, mem)
                for pn, (y0, y1) in PERIODS.items():
                    grids[f'{ssp}_{pn}'] = (clim(f, y0, y1) - h).transpose('month', 'lat', 'lon').astype('float32')
            xr.Dataset(grids).to_netcdf(gf, engine='h5netcdf')
        for si, ssp in enumerate(SSPS):
            for pi, pn in enumerate(PERIODS):
                dlt = grids[f'{ssp}_{pn}']
                a = dlt.values.reshape(12, -1)
                lon2, lat2 = np.meshgrid(dlt.lon.values, dlt.lat.values)
                ok = np.isfinite(a).all(0)
                tree = cKDTree(xyz(lon2.ravel()[ok], lat2.ravel()[ok]))
                dist, idx = tree.query(q, k=4)
                wt = 1 / np.maximum(dist, 1e-9); wt /= wt.sum(1, keepdims=True)
                res[:, gi, si, pi, :] = np.einsum('ck,mck->cm', wt, a[:, ok][:, idx])
                print(f'  {ssp} {pn}: rice-cell mean dT {np.nanmean(res[:, gi, si, pi]):.2f} K', flush=True)
    ds = xr.Dataset(
        {'dtsl': (('cell', 'gcm', 'ssp', 'period', 'month'), res,
                  {'units': 'K', 'long_name': f'CMIP6 tsl change at {ZT} m vs {HIST[0]}-{HIST[1]}'})},
        coords={'cell': cells.cell.values, 'lon': ('cell', cells.lon.values), 'lat': ('cell', cells.lat.values),
                'gcm': list(GCMS), 'ssp': SSPS, 'period': list(PERIODS), 'month': np.arange(1, 13)})
    enc = {'dtsl': {'compression': 'gzip', 'compression_opts': 4}}
    ds.to_netcdf(f'{out}/tsl_delta_cells.nc', engine='h5netcdf', encoding=enc)
    print('wrote', out)

if __name__ == '__main__':
    main(*sys.argv[1:4])
