"""Combine the ERA5-Land baseline with CMIP6 changes into one long per-cell monthly forcing table.

Baseline (gcm='ERA5-Land', ssp='historical', period='1995-2014') comes from era5land_clim_cells.csv.
Each future case (gcm x ssp x period) is baseline + CMIP6 change from the same GCM vs its own 1995-2014:
  tsoil_C  = baseline + dtsl (tsl at 0.10 m, additive)
  p_mm, et_mm, ro_mm = baseline * ratio of future/historical pr, evspsbl, mrro (ratio clipped to 0.25-4;
                       months where the GCM historical value is below 0.1 mm/day use the additive change)
  swvl1    = baseline * ratio of mrsos (clipped 0.5-2), and wsat = swvl1 / theta_sat clipped to 0-1

Usage: python3 combine_forcing.py <era5land_clim_cells.csv> <tsl_delta_cells.nc> <water_clim_cells.nc> <out.parquet>
"""
import sys
import numpy as np, pandas as pd, xarray as xr

SPD = 86400.0

def main(basef, tslf, waterf, outf):
    b = pd.read_csv(basef).sort_values(['cell', 'month']).reset_index(drop=True)
    cells = b.cell.values.reshape(-1, 12)[:, 0]
    nc = len(cells)
    B = {k: b[k].values.reshape(nc, 12) for k in ('tsoil_C', 'swvl1', 'p_mm', 'et_mm', 'ro_mm', 'theta_sat')}
    t = xr.open_dataset(tslf, engine='h5netcdf').dtsl.sel(cell=cells)
    w = xr.open_dataset(waterf, engine='h5netcdf').sel(cell=cells)
    days = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    out = [b.assign(gcm='ERA5-Land', ssp='historical', period='1995-2014')]
    for g in w.gcm.values:
        for s in w.ssp.values:
            for p in w.period.values:
                r = {'cell': np.repeat(cells, 12), 'month': np.tile(np.arange(1, 13), nc)}
                r['tsoil_C'] = (B['tsoil_C'] + t.sel(gcm=g, ssp=s, period=p).transpose('cell', 'month').values).ravel()
                for var, col in (('pr', 'p_mm'), ('evspsbl', 'et_mm'), ('mrro', 'ro_mm')):
                    h = w[f'{var}_hist'].sel(gcm=g).transpose('cell', 'month').values
                    f = w[f'{var}_fut'].sel(gcm=g, ssp=s, period=p).transpose('cell', 'month').values
                    ratio = np.clip(f / np.where(h == 0, np.nan, h), 0.25, 4)
                    add = (f - h) * SPD * days                       # mm/month
                    small = ~(np.abs(h) * SPD >= 0.1) | ~np.isfinite(ratio)
                    v = np.where(small, B[col] + add, B[col] * ratio)
                    if col != 'et_mm':
                        v = np.maximum(v, 0)
                    r[col] = v.ravel()
                h = w['mrsos_hist'].sel(gcm=g).transpose('cell', 'month').values
                f = w['mrsos_fut'].sel(gcm=g, ssp=s, period=p).transpose('cell', 'month').values
                sw = B['swvl1'] * np.nan_to_num(np.clip(f / h, 0.5, 2), nan=1.0)
                r['swvl1'] = sw.ravel()
                r['theta_sat'] = np.repeat(B['theta_sat'][:, 0], 12)
                r['wsat'] = np.clip(sw / B['theta_sat'], 0, 1).ravel()
                out.append(pd.DataFrame(r).assign(gcm=g, ssp=s, period=p))
                print(g, s, p, 'mean dT %.2f' % np.nanmean(r['tsoil_C'] - B['tsoil_C'].ravel()), flush=True)
    df = pd.concat(out, ignore_index=True)
    cols = ['cell', 'month', 'gcm', 'ssp', 'period', 'tsoil_C', 'wsat', 'p_mm', 'et_mm', 'ro_mm', 'swvl1', 'theta_sat']
    df = df[cols]
    for c in cols[5:]:
        df[c] = df[c].astype('float32')
    for c in ('gcm', 'ssp', 'period'):
        df[c] = df[c].astype('category')
    df.to_parquet(outf, index=False, compression='zstd')
    print('wrote', outf, len(df), 'rows')

if __name__ == '__main__':
    main(*sys.argv[1:5])
