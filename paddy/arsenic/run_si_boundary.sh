#!/bin/bash
# Targeted test of the predicted suppression/mobilisation boundary in baseline Si x f (paddy/docs/arsenic_step1.md,
# "Si x f boundary test"). Each case is one soil (baseline porewater Si set by irrigation-water Si and phytolith
# silica; ferrihydrite share f set by ferrihydrite vs goethite), taken through the sweep workflow (drained spin-up,
# settle without rock) and then run as a control plus a basalt dose ladder from the same settled state. Every arm uses
# the same parameters (year-matched control). Real basalt dissolution supplies the Si increment; the realised baseline
# Si, Si increment and f are read from the output, not assumed.
#
# Usage (repository root, SCEPTER_OUT set, scepter_paddy_fast built from this branch):
#   paddy/arsenic/run_si_boundary.sh <case 0-7> [settle years, default 5] [chain years, default 3]
# Output: $SCEPTER_OUT/sib_<label>/paddy, so as_sweep_summary.py runs unchanged, once per dose:
#   python3 paddy/arsenic/as_sweep_summary.py $SCEPTER_OUT/sib_<label>/paddy <dose>
#
# Fresh-season check (0.3 m, 3% ferrihydrite): porewater Si ~ irrigation Si + 10-60 uM (0.03/0.08/0.2/0.35 mM gave
# 36/97/251/396 uM); settled sweep soils sit ~12% lower (0.2 mM gave ~220 uM). Override with IRR050 etc.
# Doses 1/4/20 t/ha aim at +~125/+500/+2500 uM Si (sweep: 20 t/ha gave ~+2500 at 220 uM baseline).
# Cases: label  irrigation Si (mM)  phytolith silica (wt)  ferrihydrite (wt)  goethite (wt)
#   target baseline Si ~50 uM (f high, very high), ~100 uM (f 0.8-0.95), ~300 uM (f 0.9-0.99), ~220 uM sweep reference
CASES=(
  "si050_fhi  ${IRR050:-0.05} ${AMS050:-0.001} 0.03 0.010"
  "si050_fvh  ${IRR050:-0.05} ${AMS050:-0.001} 0.03 0.002"
  "si100_fmd  ${IRR100:-0.09} ${AMS100:-0.003} 0.03 0.015"
  "si100_fhi  ${IRR100:-0.09} ${AMS100:-0.003} 0.03 0.004"
  "si300_fhi  ${IRR300:-0.30} ${AMS300:-0.01}  0.03 0.004"
  "si300_fvh  ${IRR300:-0.30} ${AMS300:-0.01}  0.03 0.0008"
  "si220_fhi  0.2 0.01 0.03 0.003"
  "si220_fvh  0.2 0.01 0.03 0.0008"
)
DOSES=(${DOSES:-1 4 20})
set -u
i=$1; NSET=${2:-5}; NCH=${3:-3}
read -r lab irr ams f g <<< "${CASES[$i]}"
tag=as_f${f}_g${g}_a1
AS_SOIL=${AS_SOIL:-10}; AS_IRR=${AS_IRR:-100}
source stanage_handoff/v22_options.env
export XRSS=,as3up KEEPPH=1 XRS=,fhage
EXE=${AS_EXE:-scepter_paddy_fast}; K=${AS_K:-paddy/krxn_calibrated_v25_cal.in}
PROD=$(echo $PROD | sed "s#--k=[^ ]*#--k=$K#")
S="--fe3=$f --gt=$g --irralk=1 --irrsi=$irr --amsi=$ams --as=$AS_SOIL:$AS_IRR:0.1 --exe=$EXE ${SOILX:-}"
OUT=${SCEPTER_OUT%/}/sib_$lab/; mkdir -p $OUT/paddy
export SCEPTER_OUT=$OUT
ST=${ST:-status/si_boundary_$lab.txt}
grep -q "^${tag}_spin rc 0" $ST 2>/dev/null || {
  python3 paddy/run_flooded_column.py ${tag}_spin fe2o2,amo2o,as3sp,as5sp,as5rd,as3ox,as5ra,as3oa 3 --fallow --basalt=0 $SOIL $P21 $S \
    > /dev/null 2>&1; echo "${tag}_spin rc $?" >> $ST; }
paddy/multiyear/run_chain.sh ${tag}_set $NSET 0 ${tag}_spin $ST $PROD $P21 $S --kaq=0.74 || exit 1
for d in 0 ${DOSES[@]}; do
  [ -d $OUT/paddy/flood_${tag}_st$d ] || cp -a $OUT/paddy/flood_${tag}_set_y${NSET}f $OUT/paddy/flood_${tag}_st$d
  NOPRUNE=1 paddy/multiyear/run_chain.sh ${tag}_d$d $NCH $d ${tag}_st$d $ST $PROD $P21 $S --kaq=0.74 &
done; wait
