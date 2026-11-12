#!/bin/bash
# Stress test of the arsenic response map: rerun one sweep soil from its settled state (run_as_sweep.sh, settle year
# NSET) for one control + basalt year under a changed assumption. Output goes to $SCEPTER_OUT/stress_<scenario>/paddy,
# so as_sweep_summary.py / as_sweep_years.py run unchanged on that directory (pass the dose).
# Usage (repository root, SCEPTER_OUT as for the sweep): paddy/arsenic/run_as_stress.sh <soil index 0-26> <scenario> [NSET, default 5]
# Scenarios: dose5, dose10, dose40 (basalt t/ha); ox01 (As(III) oxidation as3ox/as3oa 10x slower, both arms);
#            fsite2 (twice the ferrihydrite sites available to As, both arms)
#            worst40 (ox01 + fsite2 together, 40 t/ha: combined worst case)
# Every scenario runs its own control from the same settled state with the same parameters (year-matched control).
set -u
i=$1; sc=$2; NSET=${3:-5}
FE3=(0.004 0.012 0.03); GT=(0.005 0.02 0.06); ALK=(0 1 5)
f=${FE3[$((i/9))]} g=${GT[$(((i/3)%3))]} a=${ALK[$((i%3))]}; tag=as_f${f}_g${g}_a${a}
source stanage_handoff/v22_options.env
export XRSS=,as3up KEEPPH=1 XRS=,fhage NOPRUNE=1
EXE=${AS_EXE:-scepter_paddy_fast}; K=${AS_K:-paddy/krxn_calibrated_v25_cal.in}; DOSE=20
S="--fe3=$f --gt=$g --irralk=$a --as=${AS_SOIL:-10}:${AS_IRR:-100}:0.1 --exe=$EXE ${SOILX:-}"
case $sc in
  dose5) DOSE=5;; dose10) DOSE=10;; dose40) DOSE=40;;
  ox01) K=paddy/arsenic/krxn_v25_as3ox_slow.in;;
  fsite2) S="$S --assite=2:0.02";;
  worst40) DOSE=40; K=paddy/arsenic/krxn_v25_as3ox_slow.in; S="$S --assite=2:0.02";;
  *) echo "unknown scenario $sc" >&2; exit 2;;
esac
PROD=$(echo $PROD | sed "s#--k=[^ ]*#--k=$K#")
BASE=$SCEPTER_OUT; OUT=${BASE%/}/stress_$sc/; mkdir -p $OUT/paddy
ST=${ST:-status/as_stress_${i}_$sc.txt}
src=${BASE%/}/paddy/flood_${tag}_set_y${NSET}f
[ -d $src ] || { echo "no settled state $src" >&2; exit 1; }
for d in 0 $DOSE; do
  [ -d $OUT/paddy/flood_${tag}_st$d ] || cp -a $src $OUT/paddy/flood_${tag}_st$d
  SCEPTER_OUT=$OUT paddy/multiyear/run_chain.sh ${tag}_d$d 1 $d ${tag}_st$d $ST $PROD $P21 $S --kaq=0.74 &
done; wait
