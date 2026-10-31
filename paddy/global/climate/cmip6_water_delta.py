"""Monthly CMIP6 water-variable climatologies at rice cells (companion to cmip6_tsl_delta.py).

Variables: mrsos (soil water in top 10 cm, kg m-2), pr and evspsbl (kg m-2 s-1), mrro (kg m-2 s-1).
For each GCM: historical 1995-2014 monthly mean, and for each SSP x period the future 20-yr monthly mean,
on the native grid (land points only), then the inverse-distance mean of the 4 nearest land points at
each rice cell. GFDL-ESM4 has no historical mrsos in the Pangeo archive, so its historical top-10 cm
water is mrsol summed over the layers above 0.10 m (the definition of mrsos).
Source: Pangeo CMIP6 Google Cloud archive (gs://cmip6, anonymous).

Usage: python3 cmip6_water_delta.py <pangeo-cmip6.csv> <cells.csv> <out_dir>
Writes water_clim_cells.nc with <var>_hist (cell x gcm x month) and <var>_fut (cell x gcm x ssp x period x month),
and grid_water_<gcm>.nc (the same climatologies on each GCM grid). Resumes from grid files.
"""
import sys, os
import numpy as np, pandas as pd, xarray as xr
from scipy.spatial import cKDTree

GCMS = {'GFDL-ESM4': 'r1i1p1f1', 'IPSL-CM6A-LR': 'r1i1p1f1', 'MIROC6': 'r1i1p1f1',
        'BCC-CSM2-MR': 'r1i1p1f1', 'UKESM1-0-LL': 'r1i1p1f2'}
SSPS = ['ssp245', 'ssp585']
PERIODS = {'2041-2060': (2041, 2060), '2061-2080': (2061, 2080)}
HIST = (1995, 2014)
VARS = {'mrsos': 'Lmon', 'pr': 'Amon', 'evspsbl': 'Amon', 'mrro': 'Lmon'}

def xyz(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]

def open_var(cat, var, table, gcm, exp, mem):
    s = cat[(cat.variable_id == var) & (cat.table_id == table) & (cat.source_id == gcm)
            & (cat.experiment_id == exp) & (cat.member_id == mem)]
    z = s.sort_values('version').zstore.iloc[-1]
    print(' ', z, flush=True)
    return xr.open_zarr(z, storage_options={'token': 'anon'}, consolidated=True)[var]

def sel(v, y0, y1):
    yr = v.time.dt.year
    v = v.isel(time=((yr >= y0) & (yr <= y1)).values)
    assert v.time.size == 12 * (y1 - y0 + 1), (v.time.size, y0, y1)
    return v

def clim(v, y0, y1):
    return sel(v, y0, y1).load().groupby('time.month').mean('time')

def gfdl_hist_mrsos(cat, mem):
    v = open_var(cat, 'mrsol', 'Emon', 'GFDL-ESM4', 'historical', mem)
    dz = [d for d in v.dims if d not in ('time', 'lat', 'lon')][0]
    v = v.rename({dz: 'depth'})
    d = v.depth.values
    # layer bottoms: midpoints between centres, extrapolated for the first layer
    top = np.r_[0, (d[:-1] + d[1:]) / 2]
    bot = np.r_[(d[:-1] + d[1:]) / 2, d[-1] + (d[-1] - d[-2]) / 2]
    frac = np.clip((0.10 - top) / (bot - top), 0, 1)
    print('   GFDL mrsol layers', np.round(d, 3), 'weights', np.round(frac, 2), flush=True)
    k = np.nonzero(frac)[0]
    v = sel(v, *HIST).isel(depth=k).load()
    return (v * xr.DataArray(frac[k], dims='depth')).sum('depth').groupby('time.month').mean('time')

def main(catf, cellf, out):
    os.makedirs(out, exist_ok=True)
    cat = pd.read_csv(catf)
    cells = pd.read_csv(cellf, usecols=['cell', 'lon', 'lat'])
    q = xyz(cells.lon.values, cells.lat.values)
    nc, ng, ns, npd = len(cells), len(GCMS), len(SSPS), len(PERIODS)
    hist = {v: np.full((nc, ng, 12), np.nan, np.float32) for v in VARS}
    fut = {v: np.full((nc, ng, ns, npd, 12), np.nan, np.float32) for v in VARS}
    for gi, (gcm, mem) in enumerate(GCMS.items()):
        print(gcm, flush=True)
        gf = f'{out}/grid_water_{gcm}.nc'
        if os.path.exists(gf):
            grids = xr.open_dataset(gf, engine='h5netcdf').load()
        else:
            g = {}
            for var, tab in VARS.items():
                if var == 'mrsos' and gcm == 'GFDL-ESM4':
                    g[f'{var}_hist'] = gfdl_hist_mrsos(cat, mem)
                else:
                    g[f'{var}_hist'] = clim(open_var(cat, var, tab, gcm, 'historical', mem), *HIST)
                for ssp in SSPS:
                    f = open_var(cat, var, tab, gcm, ssp, mem)
                    for pn, (y0, y1) in PERIODS.items():
                        g[f'{var}_{ssp}_{pn}'] = clim(f, y0, y1)
            # atmosphere and land tables can carry float-noise differences in lat/lon (e.g. BCC-CSM2-MR);
            # snap every field onto the grid of the first one so the Dataset does not outer-join them
            ref = g['mrsos_hist']
            g = {k: x.reindex(lat=ref.lat, lon=ref.lon, method='nearest', tolerance=1e-3) for k, x in g.items()}
            grids = xr.Dataset({k: x.transpose('month', 'lat', 'lon').astype('float32').drop_vars(
                [c for c in x.coords if c not in ('month', 'lat', 'lon')]) for k, x in g.items()})
            grids.to_netcdf(gf, engine='h5netcdf')
        lon2, lat2 = np.meshgrid(grids.lon.values, grids.lat.values)
        # land mask: points where soil moisture is defined in every field (pr/evspsbl are global)
        ok = np.ones(lon2.size, bool)
        for k in grids.data_vars:
            if k.startswith(('mrsos', 'mrro')):
                ok &= np.isfinite(grids[k].values.reshape(12, -1)).all(0)
        tree = cKDTree(xyz(lon2.ravel()[ok], lat2.ravel()[ok]))
        dist, idx = tree.query(q, k=4)
        wt = 1 / np.maximum(dist, 1e-9); wt /= wt.sum(1, keepdims=True)
        interp = lambda k: np.einsum('ck,mck->cm', wt, grids[k].values.reshape(12, -1)[:, ok][:, idx])
        for var in VARS:
            hist[var][:, gi] = interp(f'{var}_hist')
            for si, ssp in enumerate(SSPS):
                for pi, pn in enumerate(PERIODS):
                    fut[var][:, gi, si, pi] = interp(f'{var}_{ssp}_{pn}')
            h = np.nanmean(hist[var][:, gi]); f = np.nanmean(fut[var][:, gi, 1, 1])
            print(f'  {var}: hist {h:.4g}  ssp585 2061-80 {f:.4g}  ratio {f / h:.3f}', flush=True)
    units = {'mrsos': 'kg m-2 (top 0.1 m)', 'pr': 'kg m-2 s-1', 'evspsbl': 'kg m-2 s-1', 'mrro': 'kg m-2 s-1'}
    dv = {}
    for v in VARS:
        dv[f'{v}_hist'] = (('cell', 'gcm', 'month'), hist[v], {'units': units[v], 'period': f'{HIST[0]}-{HIST[1]}'})
        dv[f'{v}_fut'] = (('cell', 'gcm', 'ssp', 'period', 'month'), fut[v], {'units': units[v]})
    ds = xr.Dataset(dv, coords={'cell': cells.cell.values, 'lon': ('cell', cells.lon.values),
                                'lat': ('cell', cells.lat.values), 'gcm': list(GCMS), 'ssp': SSPS,
                                'period': list(PERIODS), 'month': np.arange(1, 13)})
    enc = {k: {'compression': 'gzip', 'compression_opts': 4} for k in dv}
    ds.to_netcdf(f'{out}/water_clim_cells.nc', engine='h5netcdf', encoding=enc)
    print('wrote', out)

if __name__ == '__main__':
    main(*sys.argv[1:4])
