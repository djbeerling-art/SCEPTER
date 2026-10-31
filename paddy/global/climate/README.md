# Monthly soil climate forcing for the riceland runs

Target (agreed with the Stanage thread): one long table keyed on the `cell` id of `data/global/cells.csv`
(branch claude/project-thread-lst1mi, 109,319 rice cells within 250 km of basalt), columns
`cell, month, gcm, ssp, period, tsoil_C, wsat, p_mm, et_mm` (+ `ro_mm, swvl1, theta_sat`), for the
baseline and each CMIP6 case. The Stanage thread turns it into T_temp.in etc. via `--clim`.

## Status (2026-10-31)

| Piece | State | Files |
|---|---|---|
| CMIP6 soil-T change (tsl, 0.10 m) | done | `../soil-temperature/cmip6_tsl_delta/tsl_delta_cells.nc` |
| CMIP6 water climatologies (mrsos, pr, evspsbl, mrro), hist 1995-2014 and futures | done | `cmip6_water/water_clim_cells.nc`, `cmip6_water/grid_water_<GCM>.nc` |
| WorldClim 2.1 5' tmin/tmax/prec at cells (1970-2000 + CMIP6) | done | `worldclim_cells.parquet` |
| ERA5-Land 1995-2014 baseline (stl1/2, swvl1/2, tp, e, ro, t2m) | **waiting for CDS key** | `era5land_baseline.py` |
| Final table | script ready, needs baseline | `combine_forcing.py` -> `forcing_cells.parquet` |

## Method
- GCMs (r1, f2 for UKESM): BCC-CSM2-MR, GFDL-ESM4, IPSL-CM6A-LR, MIROC6, UKESM1-0-LL; SSP2-4.5 and SSP5-8.5;
  2041-2060 and 2061-2080, each vs the same model's 1995-2014. Pangeo gs://cmip6. Rice cells take the
  inverse-distance mean of the 4 nearest GCM land points.
- GFDL-ESM4 has no historical `mrsos` on Pangeo; its historical top-10 cm water is `mrsol` summed over
  the layers above 0.10 m (future `mrsos` vs this agree to ~2%).
- Baseline: ERA5-Land monthly means, 0.1 deg, nearest land point. `tsoil_C` = stl1 (0-7 cm).
  `wsat` = swvl1 / theta_sat with theta_sat from the HTESSEL soil type (ERA5 `soil_type`, `era5_soil_type.nc`).
  `p_mm`, `et_mm` (= -e), `ro_mm` are monthly totals. Note ERA5-Land does not know about paddy flooding,
  so wsat is the rain-fed soil state; irrigated/flooded periods should be set by the crop calendar.
- Futures: tsoil + dtsl; p/et/ro × CMIP6 future/hist ratio (clipped 0.25-4, additive where hist < 0.1 mm/day);
  swvl1 × mrsos ratio (clipped 0.5-2).

## Rice-cell mean CMIP6 change, SSP5-8.5 2061-2080 (ratio future/hist)

| GCM | mrsos | pr | evspsbl | mrro |
|---|---|---|---|---|
| GFDL-ESM4 | 1.02 | 1.00 | 0.98 | 1.07 |
| IPSL-CM6A-LR | 0.98 | 1.05 | 1.00 | 1.15 |
| MIROC6 | 0.98 | 1.04 | 1.03 | 1.06 |
| BCC-CSM2-MR | 0.99 | 1.03 | 1.02 | 1.06 |
| UKESM1-0-LL | 0.92 | 1.06 | 1.01 | 1.19 |

## WorldClim notes
WorldClim's server has no GFDL-ESM4 SSP2-4.5 and only precipitation for GFDL-ESM4 SSP5-8.5, so GFDL is
absent from `worldclim_cells.parquet`; the CMIP6 fields above cover it. 1,298 cells (coastal) are NaN in
the WorldClim future grids and 12 in the baseline. BCC-CSM2-MR has no bioclim files.

## Running the baseline
```
pip install cdsapi xarray netCDF4 scipy pandas pyarrow h5netcdf
export CDSAPI_KEY=<CDS personal access token>   # and accept the ERA5-Land licence on the CDS site
python3 era5land_baseline.py cells.csv era5_soil_type.nc era5land_out 4
python3 combine_forcing.py era5land_out/era5land_clim_cells.csv \
    ../soil-temperature/cmip6_tsl_delta/tsl_delta_cells.nc cmip6_water/water_clim_cells.nc forcing_cells.parquet
```

## Running the baseline on Stanage (the CDS key is in ~/.cdsapirc there)
```
module load Python/3.11.3-GCCcore-12.3.0
pip install --user cdsapi xarray netCDF4 scipy pandas
cd <repo> && python3 paddy/global/climate/era5land_baseline.py data/global/cells.csv \
    paddy/global/climate/era5_soil_type.nc ~/fastdata/riceland/inputs/era5land 4
```
Output `era5land_clim_cells.csv.gz` (~1.3M rows, tens of MB) can be pushed back to GitHub for the combine step.
Downloads are one year (12 months x 8 variables over 41S-53N) per CDS request; raw files are deleted after
sampling, so disk use stays at a few GB.
