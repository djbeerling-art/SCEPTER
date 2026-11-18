#!/bin/bash
# Arsenic response sweep (Step 1): soil-dependent response of porewater As(III), Si and As(III) uptake to basalt.
# For each soil (reactive Fe x crystalline Fe x irrigation alkalinity) it runs a drained spin-up, a settle of season +
# fallow years without rock, then control and basalt chains from the same settled state (as the v22/v25 workflow in
# paddy/docs/USER_MANUAL.md 5.1). Arsenic is tracked from the spin-up on, so restarts stay compatible.
#
# Usage (repository root, SCEPTER_OUT set, scepter_paddy_fast built from makefile_paddy_fast):
#   paddy/arsenic/run_as_sweep.sh <soil index 0-26 | all> [settle years, default 5] [chain years, default 3] [dose t/ha, default 20]
# Stanage array:  #SBATCH --array=0-26  ...  paddy/arsenic/run_as_sweep.sh $SLURM_ARRAY_TASK_ID
# Then:           python3 paddy/arsenic/as_sweep_summary.py $SCEPTER_OUT/paddy > as_sweep_summary.csv
#
# Soil axes (edit here):
#   FE3   ferrihydrite (amfe3, wt fraction; ~0.4, 1.2, 3% = 45, 135, 340 mmol Fe/kg; oxalate Fe in S&F 2012: 54-537)
#   GT    goethite (wt fraction), the Si-insensitive As pool (0.02 sites per mol Fe; DCB - oxalate Fe)
#   ALK   --irralk (irrigation alkalinity factor; 0 = acid, 1 = reference, 5 = alkaline groundwater) -> initial pH
# SOILX (environment) appends soil options to every run, e.g. "--soilcc=0.02" for a calcareous soil.
# Arsenic (fixed): labile soil As AS_SOIL mg/kg, irrigation As AS_IRR ug/L (Bengal-type groundwater), 90% As(III).
set -u
idx=${1:-all}; NSET=${2:-5}; NCH=${3:-3}; DOSE=${4:-20}
FE3=(0.004 0.012 0.03); GT=(0.005 0.02 0.06); ALK=(0 1 5)
AS_SOIL=${AS_SOIL:-10}; AS_IRR=${AS_IRR:-100}
ST=${ST:-status_as_sweep.txt}
source stanage_handoff/v22_options.env
export XRSS=,as3up KEEPPH=1 XRS=,fhage   # v25: ferrihydrite ageing on; As(III) uptake in the season only
EXE=${AS_EXE:-scepter_paddy_fast}         # must be built from this (arsenic) branch
K=${AS_K:-paddy/krxn_calibrated_v25_cal.in}
PROD=$(echo $PROD | sed "s#--k=[^ ]*#--k=$K#")   # one --k only (run_chain.sh passes it to the fallow)

run_soil(){
  local i=$1
  local f=${FE3[$((i/9))]} g=${GT[$(((i/3)%3))]} a=${ALK[$((i%3))]}
  local tag=as_f${f}_g${g}_a${a}
  local S="--fe3=$f --gt=$g --irralk=$a --as=$AS_SOIL:$AS_IRR:0.1 --exe=$EXE ${SOILX:-}"
  # 1. drained spin-up, 3 yr
  grep -q "^${tag}_spin rc 0" $ST 2>/dev/null || {
    python3 paddy/run_flooded_column.py ${tag}_spin fe2o2,amo2o,as3sp,as5sp,as5rd,as3ox,as5ra,as3oa 3 --fallow --basalt=0 $SOIL $P21 $S \
      > /dev/null 2>&1; echo "${tag}_spin rc $?" >> $ST; }
  # 2. settle without rock
  paddy/multiyear/run_chain.sh ${tag}_set $NSET 0 ${tag}_spin $ST $PROD $P21 $S --kaq=0.74 || return 1
  # 3. control and basalt from the settled state (separate copies: tillage acts on each)
  for d in 0 $DOSE; do
    [ -d $SCEPTER_OUT/paddy/flood_${tag}_st$d ] || cp -a $SCEPTER_OUT/paddy/flood_${tag}_set_y${NSET}f $SCEPTER_OUT/paddy/flood_${tag}_st$d
    NOPRUNE=1 paddy/multiyear/run_chain.sh ${tag}_d$d $NCH $d ${tag}_st$d $ST $PROD $P21 $S --kaq=0.74 &
  done; wait
}

if [ "$idx" = all ]; then for i in $(seq 0 26); do run_soil $i; done; else run_soil $idx; fi
