"""ERA5-Land monthly 1995-2014 baseline climatology at the rice cells.

Downloads ERA5-Land monthly means (CDS dataset reanalysis-era5-land-monthly-means) one year per request,
samples each month at the rice cells, accumulates the 20-yr monthly climatology and deletes the raw file.
Each cell takes the nearest ERA5-Land land point (0.1 deg; searches outward if the nearest is sea).

Variables: stl1 (0-7 cm soil T), stl2 (7-28 cm), swvl1/swvl2 (volumetric water, m3 m-3), tp, e, ro
(monthly-mean daily totals, m/day; e negative = evaporation), t2m.
Porosity: HTESSEL saturation by ERA5 soil type (era5_soil_type.nc, 0.25 deg, from the Google ARCO ERA5
archive), used for wsat = swvl1 / theta_sat.

Needs CDSAPI_KEY (CDS personal access token) or ~/.cdsapirc, and the ERA5-Land licence accepted on the CDS website.
Usage: python3 era5land_baseline.py <cells.csv (data/global/cells.csv)> <era5_soil_type.nc> <out_dir> [n_parallel=4]
Writes era5land_clim_cells.csv.gz (cell, month, ...) and keeps per-year partial sums so it can resume.
"""
import sys, os, calendar, threading, concurrent.futures as cf
import numpy as np, pandas as pd, xarray as xr
from scipy.spatial import cKDTree

YEARS = range(1995, 2015)
VARS = {'soil_temperature_level_1': 'stl1', 'soil_temperature_level_2': 'stl2',
        'volumetric_soil_water_layer_1': 'swvl1', 'volumetric_soil_water_layer_2': 'swvl2',
        'total_precipitation': 'tp', 'total_evaporation': 'e', 'runoff': 'ro', '2m_temperature': 't2m'}
# HTESSEL saturated volumetric water by soil type 1-7 (coarse, medium, medium-fine, fine, very fine,
# organic, tropical organic); ECMWF IFS documentation Part IV, Table 8.9
THETA_SAT = {1: 0.403, 2: 0.439, 3: 0.430, 4: 0.520, 5: 0.614, 6: 0.766, 7: 0.472}
AREA = [53, -125, -41, 180]  # N, W, S, E around all rice cells
lock = threading.Lock()

def fetch(year, path):
    import cdsapi
    if os.environ.get('CDSAPI_KEY'):
        c = cdsapi.Client(url='https://cds.climate.copernicus.eu/api', key=os.environ['CDSAPI_KEY'], quiet=True)
    else:
        c = cdsapi.Client(quiet=True)   # reads ~/.cdsapirc
    c.retrieve('reanalysis-era5-land-monthly-means', {
        'product_type': ['monthly_averaged_reanalysis'], 'variable': list(VARS),
        'year': [str(year)], 'month': [f'{m:02d}' for m in range(1, 13)], 'time': ['00:00'],
        'area': AREA, 'data_format': 'netcdf', 'download_format': 'unarchived'}, path)

def nearest_land(lon_g, lat_g, valid, lon, lat):
    """Index (into the flattened grid) of the nearest valid point for each cell."""
    LO, LA = np.meshgrid(lon_g, lat_g)
    ok = np.flatnonzero(valid.ravel())
    lo = np.where(LO.ravel()[ok] > 180, LO.ravel()[ok] - 360, LO.ravel()[ok])
    xyz = lambda a, b: np.c_[np.cos(np.radians(b)) * np.cos(np.radians(a)),
                             np.cos(np.radians(b)) * np.sin(np.radians(a)), np.sin(np.radians(b))]
    d, i = cKDTree(xyz(lo, LA.ravel()[ok])).query(xyz(lon, lat))
    return ok[i], d * 6371.0

def main(cellf, sltf, out, npar=4):
    os.makedirs(f'{out}/years', exist_ok=True)
    cells = pd.read_csv(cellf, usecols=['cell', 'lon', 'lat'])
    state = {}

    def one(year):
        part = f'{out}/years/{year}.npz'
        if os.path.exists(part):
            return year
        raw = f'{out}/raw_{year}.nc'
        if not os.path.exists(raw):
            fetch(year, raw + '.tmp'); os.replace(raw + '.tmp', raw)
        ds = xr.open_dataset(raw)
        tdim = 'valid_time' if 'valid_time' in ds.dims else 'time'
        ds = ds.sortby(tdim)
        assert ds[tdim].size == 12, ds[tdim].size
        with lock:
            if 'idx' not in state:
                valid = np.isfinite(ds['swvl1'].isel({tdim: 6}).values) & np.isfinite(ds['stl1'].isel({tdim: 6}).values)
                state['idx'], dist = nearest_land(ds.longitude.values, ds.latitude.values, valid,
                                                  cells.lon.values, cells.lat.values)
                state['shape'] = valid.shape
                np.save(f'{out}/era5land_index.npy', state['idx'])
                print(f'grid {valid.shape}; cell->land point distance: median {np.median(dist):.1f} km, '
                      f'max {dist.max():.1f} km, >10 km: {(dist > 10).sum()}', flush=True)
        idx = state['idx']
        arr = {v: ds[v].values.reshape(12, -1)[:, idx].astype('float32') for v in VARS.values()}
        np.savez(part, **arr)
        ds.close(); os.remove(raw)
        print(year, 'done', flush=True)
        return year

    with cf.ThreadPoolExecutor(int(npar)) as ex:
        for f in [ex.submit(one, y) for y in YEARS]:
            f.result()

    acc = {v: np.zeros((12, len(cells)), 'float64') for v in VARS.values()}
    days = np.zeros(12)
    for y in YEARS:
        z = np.load(f'{out}/years/{y}.npz')
        dim = np.array([calendar.monthrange(y, m)[1] for m in range(1, 13)])
        for v in acc:
            x = z[v].astype('float64')
            if v in ('tp', 'e', 'ro'):
                x = x * 1000 * dim[:, None]     # m/day monthly mean -> mm per month
            acc[v] += x
    n = len(YEARS)
    clim = {v: (a / n).astype('float32') for v, a in acc.items()}

    st = xr.open_dataset(sltf).soil_type
    st = st.assign_coords(longitude=((st.longitude + 180) % 360) - 180).sortby('longitude')
    sv = np.round(st.values)
    j, _ = nearest_land(st.longitude.values, st.latitude.values, sv > 0, cells.lon.values, cells.lat.values)
    slt = sv.ravel()[j].astype(int)
    tsat = np.vectorize(THETA_SAT.get)(slt).astype('float32')

    rows = pd.DataFrame({
        'cell': np.tile(cells.cell.values, 12), 'month': np.repeat(np.arange(1, 13), len(cells)),
        'tsoil_C': (clim['stl1'] - 273.15).ravel(), 'tsoil2_C': (clim['stl2'] - 273.15).ravel(),
        'tair_C': (clim['t2m'] - 273.15).ravel(),
        'swvl1': clim['swvl1'].ravel(), 'swvl2': clim['swvl2'].ravel(),
        'theta_sat': np.tile(tsat, 12), 'soil_type': np.tile(slt, 12),
        'p_mm': clim['tp'].ravel(), 'et_mm': (-clim['e']).ravel(), 'ro_mm': clim['ro'].ravel()})
    rows['wsat'] = np.clip(rows.swvl1 / rows.theta_sat, 0, 1)
    rows = rows.sort_values(['cell', 'month'])
    rows.to_csv(f'{out}/era5land_clim_cells.csv.gz', index=False, float_format='%.4g')
    print(rows.describe().T.to_string(), flush=True)

if __name__ == '__main__':
    main(*sys.argv[1:])
