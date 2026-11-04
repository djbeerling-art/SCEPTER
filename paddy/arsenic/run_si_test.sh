#!/bin/bash
# Si-only arsenic test from settled flooded soils: for a sweep soil (run_as_sweep.sh index), spin up and settle, then run
# one season + fallow from the settled state with and without silica gel added to the plough layer (no basalt, so Si
# rises without the basalt pH or cation effects). Compare with paddy/arsenic/gt_si_partition.py-style diagnostics.
# Usage (repository root, SCEPTER_OUT set): paddy/arsenic/run_si_test.sh <soil index 0-26> [settle years, default 2] [silica wt, default 0.077]
set -u
i=$1; NSET=${2:-2}; SIWT=${3:-0.077}
FE3=(0.004 0.012 0.03); GT=(0.005 0.02 0.06); ALK=(0 1 5)
f=${FE3[$((i/9))]} g=${GT[$(((i/3)%3))]} a=${ALK[$((i%3))]}; tag=as_f${f}_g${g}_a${a}
ST=${ST:-status_si_test.txt}
paddy/arsenic/run_as_sweep.sh $i $NSET 0 || exit 1          # spin-up + settle only (no chains)
source stanage_handoff/v22_options.env
export XRSS=,as3up KEEPPH=1 XRS=,fhage NOPRUNE=1
EXE=${AS_EXE:-scepter_paddy_fast}; K=${AS_K:-paddy/krxn_calibrated_v25_cal.in}
PROD=$(echo $PROD | sed "s#--k=[^ ]*#--k=$K#")
S="--fe3=$f --gt=$g --irralk=$a --as=${AS_SOIL:-10}:${AS_IRR:-100}:0.1 --exe=$EXE ${SOILX:-}"
for arm in sc ss; do
  d=$SCEPTER_OUT/paddy/flood_${tag}_$arm
  [ -d $d ] || { cp -a $SCEPTER_OUT/paddy/flood_${tag}_set_y${NSET}f $d
                 [ $arm = ss ] && python3 paddy/arsenic/add_silica.py $d $SIWT ${TILL:-0.2}; }
  paddy/multiyear/run_chain.sh ${tag}_$arm 1 0 ${tag}_$arm $ST $PROD $P21 $S --kaq=0.74 &
done; wait
