#!/bin/bash
# Run the ERA5-Land 1995-2014 monthly baseline (era5land_baseline.py, from the "Build monthly soil climate forcing"
# thread) on Stanage. The CDS token is read from ~/.cdsapirc (chmod 600) into CDSAPI_KEY; it is never printed.
# Run from the repo root on the LOGIN node (downloads need internet), e.g. inside tmux or with nohup:
#   nohup bash paddy/global/forcing/run_era5land_stanage.sh > logs/era5land.log 2>&1 &
# Resumable: finished years are kept in $OUT/years/*.npz.
set -eu
source stanage_handoff/slurm/modules.sh
python3 -m pip install --user -q cdsapi xarray netCDF4 h5netcdf scipy pandas
export CDSAPI_KEY=$(awk '/^key:/{print $2}' ~/.cdsapirc)
OUT=${OUT:-/mnt/parscratch/users/$USER/riceland_inputs/era5land}
mkdir -p $OUT logs
python3 paddy/global/forcing/era5land_baseline.py data/global/cells.csv paddy/global/forcing/era5_soil_type.nc $OUT ${NPAR:-4}
ls -l $OUT/era5land_clim_cells.csv*
