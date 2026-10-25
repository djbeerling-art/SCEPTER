# Build cells.csv (the archetypes.py input) from global datasets, one row per 5-arcmin MapSPAM rice cell.
#   rice area   MapSPAM 2010 v2.0 harvested area (spam2010V2r0_global_H_TA/_TI.csv): area_irr = TI, area_rf = TA - TI
#   soil        SoilGrids 2.0, 5 km aggregates (<var>_<depth>_mean_5000.tif, Goode homolosine): phh2o, cec, soc, clay,
#               0-15 cm as the thickness-weighted mean of 0-5 and 5-15 cm; gaps filled from the nearest valid pixel
#               within FILL pixels (SoilGrids masks water and built-up land, common in deltas)
#   temp        WorldClim 2.1 5-arcmin monthly tavg: mean of the warmest NSEAS consecutive months, a stand-in for the
#               rice season where no crop calendar is used (tropics barely differ; temperate rice is a summer crop)
#   basalt      GLiM basic volcanics (xx = 'vb') as a raster or point set (--glim): cells further than DMAX km
#               from the nearest outcrop are dropped; dist_basalt (km) is written either way
#   rock        --rocks: CSV of rock recipes with lon, lat, rock (rock_recipes.py); each cell takes the nearest
#               sample's recipe. Without it every cell gets --rock=<name> (default hillhouse).
# Usage: python3 prep_cells.py <data dir> <out cells.csv> [--glim=<vb points csv: lon,lat>] [--rocks=<csv>]
#        [--rock=hillhouse] [--dmax=250] [--nseas=4] [--fill=3] [--minarea=1]
import sys, os
import numpy as np, pandas as pd, rasterio
from pyproj import Transformer
from sklearn.neighbors import BallTree

R_EARTH = 6371.

def opt(name, default):
    return next((a.split('=', 1)[1] for a in sys.argv if a.startswith(f'--{name}=')), default)

def rice_cells(d, minarea):
    cols = ['cell5m', 'iso3', 'x', 'y', 'rice_a']
    ta = pd.read_csv(f'{d}/spam/spam2010V2r0_global_H_TA.csv', usecols=cols, encoding='latin-1')
    ti = pd.read_csv(f'{d}/spam/spam2010V2r0_global_H_TI.csv', usecols=['cell5m', 'rice_i'], encoding='latin-1')
    c = ta.merge(ti, on='cell5m', how='left').fillna({'rice_i': 0.})
    c = c[c.rice_a >= minarea]
    return pd.DataFrame({'cell': c.cell5m.values, 'lon': c.x.values, 'lat': c.y.values, 'iso3': c.iso3.values,
                         'area_irr': c.rice_i.values, 'area_rf': np.maximum(c.rice_a.values - c.rice_i.values, 0.)})

def sample(path, lon, lat, fill, to_crs=None):
    """Raster value at (lon, lat); nodata filled from the nearest valid pixel within `fill` pixels (else NaN)."""
    with rasterio.open(path) as r:
        a = r.read(1).astype(float); nd = r.nodata
        if nd is not None: a[(a == nd) | (a < -1e30)] = np.nan
        x, y = (Transformer.from_crs('EPSG:4326', r.crs, always_xy=True).transform(lon, lat)
                if r.crs.to_epsg() != 4326 else (lon, lat))
        row, col = rasterio.transform.rowcol(r.transform, x, y)
    row = np.clip(np.asarray(row), 0, a.shape[0]-1); col = np.clip(np.asarray(col), 0, a.shape[1]-1)
    v = a[row, col]
    for k in range(1, fill+1):                     # rings of growing radius, nearest-ring mean
        m = np.isnan(v)
        if not m.any(): break
        acc = np.zeros(m.sum()); n = np.zeros(m.sum())
        for di in range(-k, k+1):
            for dj in range(-k, k+1):
                if max(abs(di), abs(dj)) != k: continue
                w = a[np.clip(row[m]+di, 0, a.shape[0]-1), np.clip(col[m]+dj, 0, a.shape[1]-1)]
                ok = ~np.isnan(w); acc[ok] += w[ok]; n[ok] += 1
        vv = v[m]; vv[n > 0] = acc[n > 0]/n[n > 0]; v[m] = vv
    return v

def soil(d, lon, lat, fill):
    # SoilGrids units: phh2o pH*10, cec mmol(c)/kg, soc dg/kg, clay g/kg -> pH, cmol(c)/kg, wt %, wt %
    scale = {'phh2o': 0.1, 'cec': 0.1, 'soc': 0.01, 'clay': 0.1}
    out = {}
    for v, s in scale.items():
        a, b = (sample(f'{d}/soil/{v}_{z}_mean_5000.tif', lon, lat, fill) for z in ('0-5cm', '5-15cm'))
        out[v] = s*(5*a + 10*b)/15
    return out['phh2o'], out['cec'], out['soc'], out['clay']

def season_temp(d, lon, lat, nseas, fill):
    t = np.array([sample(f'{d}/wc/wc2.1_5m_tavg_{m:02d}.tif', lon, lat, fill) for m in range(1, 13)])
    roll = np.array([np.mean(np.take(t, range(s, s+nseas), axis=0, mode='wrap'), axis=0) for s in range(12)])
    return roll.max(axis=0)

def nearest(lon, lat, plon, plat):
    """Great-circle distance (km) and index of the nearest point in (plon, plat) for every (lon, lat)."""
    tree = BallTree(np.radians(np.c_[plat, plon]), metric='haversine')
    dist, idx = tree.query(np.radians(np.c_[lat, lon]), k=1)
    return dist[:, 0]*R_EARTH, idx[:, 0]

if __name__ == '__main__':
    d, out = sys.argv[1], sys.argv[2]
    dmax, nseas, fill = float(opt('dmax', 250)), int(opt('nseas', 4)), int(opt('fill', 3))
    c = rice_cells(d, float(opt('minarea', 1)))
    print(f'{len(c)} rice cells, {c.area_irr.sum()/1e6:.1f} Mha irrigated + {c.area_rf.sum()/1e6:.1f} Mha rainfed')
    c['ph'], c['cec'], c['soc'], c['clay'] = soil(d, c.lon.values, c.lat.values, fill)
    c['temp'] = season_temp(d, c.lon.values, c.lat.values, nseas, fill)
    glim = opt('glim', None)
    if glim:
        vb = pd.read_csv(glim)
        c['dist_basalt'], _ = nearest(c.lon.values, c.lat.values, vb.lon.values, vb.lat.values)
        c = c[c.dist_basalt <= dmax]
        print(f'{len(c)} cells within {dmax:.0f} km of basic volcanics, {(c.area_irr+c.area_rf).sum()/1e6:.1f} Mha')
    rocks = opt('rocks', None)
    if rocks:
        r = pd.read_csv(rocks)
        c['dist_rock'], i = nearest(c.lon.values, c.lat.values, r.lon.values, r.lat.values)
        c['rock'] = r.rock.values[i]
    else:
        c['rock'] = opt('rock', 'hillhouse')
    bad = c[['ph', 'cec', 'soc', 'clay', 'temp']].isna().any(axis=1)
    if bad.any():
        print(f'dropping {bad.sum()} cells without soil/climate data ({(c.area_irr+c.area_rf)[bad].sum()/1e3:.0f} kha)')
        c = c[~bad]
    c.round({'lon': 4, 'lat': 4, 'area_irr': 1, 'area_rf': 1, 'ph': 2, 'cec': 1, 'soc': 3, 'clay': 1, 'temp': 2,
             'dist_basalt': 1, 'dist_rock': 1}).to_csv(out, index=False)
    print(f'wrote {out}: {len(c)} cells, {(c.area_irr+c.area_rf).sum()/1e6:.1f} Mha')
