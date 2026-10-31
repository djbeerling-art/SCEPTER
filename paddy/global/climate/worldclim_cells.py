"""Sample WorldClim 2.1 5-arcmin monthly tmin/tmax/prec (1970-2000 baseline and CMIP6 futures) at rice cells.

Usage: python3 worldclim_cells.py <cells.csv> <base_dir> <cmip6_dir> <out.parquet>
base_dir holds wc2.1_5m_<var>_<MM>.tif; cmip6_dir is the output of download_worldclim_cmip6.sh.
Output (long): cell, month, gcm, ssp, period, tmin_C, tmax_C, tavg_C, prec_mm. Baseline rows have
gcm='WorldClim', ssp='historical', period='1970-2000'. Missing future files are skipped.
"""
import sys, glob, os
import numpy as np, pandas as pd, rasterio

def sample(path, xy, band=None):
    with rasterio.open(path) as r:
        bands = [band] if band else list(range(1, r.count + 1))
        v = np.array(list(r.sample(xy, indexes=bands)), dtype='float32')
        nd = r.nodata
    if nd is not None:
        v[v == nd] = np.nan
    v[v < -1e30] = np.nan
    return v  # (ncell, nband)

def main(cellf, based, cmipd, outf):
    c = pd.read_csv(cellf, usecols=['cell', 'lon', 'lat'])
    xy = list(zip(c.lon.values, c.lat.values))
    n = len(c)
    frames = []
    base = {v: np.column_stack([sample(f'{based}/wc2.1_5m_{v}_{m:02d}.tif', xy)[:, 0] for m in range(1, 13)])
            for v in ('tmin', 'tmax', 'prec')}
    frames.append(('WorldClim', 'historical', '1970-2000', base))
    cases = {}
    for f in sorted(glob.glob(f'{cmipd}/*/*/wc2.1_5m_*.tif')):
        _, _, v, g, s, p = os.path.basename(f)[:-4].split('_')
        if v in ('tmin', 'tmax', 'prec'):
            cases.setdefault((g, s, p), {})[v] = sample(f, xy)
    for (g, s, p), d in cases.items():
        if len(d) < 3:
            print('incomplete, skipped:', g, s, p, sorted(d)); continue
        frames.append((g, s, p, d))
    out = []
    for g, s, p, d in frames:
        out.append(pd.DataFrame({'cell': np.repeat(c.cell.values, 12), 'month': np.tile(np.arange(1, 13), n),
                                 'gcm': g, 'ssp': s, 'period': p,
                                 'tmin_C': d['tmin'].ravel(), 'tmax_C': d['tmax'].ravel(),
                                 'tavg_C': ((d['tmin'] + d['tmax']) / 2).ravel(), 'prec_mm': d['prec'].ravel()}))
        print(g, s, p, 'NaN cells:', int(np.isnan(d['tmin'][:, 0]).sum()), flush=True)
    df = pd.concat(out, ignore_index=True)
    for k in ('gcm', 'ssp', 'period'):
        df[k] = df[k].astype('category')
    df.to_parquet(outf, index=False, compression='zstd')
    print('wrote', outf, len(df))

if __name__ == '__main__':
    main(*sys.argv[1:5])
