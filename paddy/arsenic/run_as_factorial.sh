#!/bin/bash
# Representative subset of the baseline Si x basalt dose x ferrihydrite share x phosphate factorial (David, 3 Nov 2026),
# chosen to validate the correction-A predictor, especially at high baseline Si and low dose. Each case is one soil
# taken through the sweep workflow (drained spin-up, settle without rock), then a control and 5/10/20 t/ha basalt
# from the same settled state, all with the same parameters (year-matched control). Realised baseline Si, Si
# increment, f and pH are read from the output (as_sweep_summary.py), not assumed.
#
# Usage (repository root, SCEPTER_OUT set):  paddy/arsenic/run_as_factorial.sh <case 0-13> [settle yrs 5] [chain yrs 3]
# Needs an executable built from this branch at or after the as_kd.in commit (AS_EXE, default scepter_paddy_fast_p),
# because phosphate enters through --askd. All 14 cases (both P levels) use that same executable.
# Output: $SCEPTER_OUT/fac_<label>/paddy; then per dose: python3 paddy/arsenic/as_sweep_summary.py <dir> <dose>
#
# Levels
#   baseline Si 100/300/700/1000 uM: irrigation Si 0.10/0.30/0.70/1.00 mM (porewater sits near irrigation Si;
#     Bangladesh shallow tubewells: median 733 uM, p10-p90 484-1043)
#   f 0.3/0.7/0.95: goethite set from the as3sp terms at that Si with 1.2% ferrihydrite (fhage lowers f a little)
#   phosphate: labile P 0.02 (low) or 0.1 (high) mol per mol ferrihydrite Fe, as refits of the DLM sorption
#     constants (kd_constants_Pr0.02/0.1.json). P is a fixed background competitor (SCEPTER carries no P), and
#     the high-P fit covers 492 of the 882 grid points (high-P, high-pH points fail to converge).
# Cases: label  irrSi(mM)  phytolith(wt)  ferrihydrite(wt)  goethite(wt)  P
CASES=(
  "s0100_f30_pL 0.10 0.003 0.012 0.0757 r0.02"
  "s0100_f95_pL 0.10 0.003 0.012 0.0017 r0.02"
  "s0300_f30_pL 0.30 0.01  0.012 0.0367 r0.02"
  "s0300_f70_pL 0.30 0.01  0.012 0.0067 r0.02"
  "s0300_f95_pL 0.30 0.01  0.012 0.0008 r0.02"
  "s0700_f30_pL 0.70 0.01  0.012 0.0192 r0.02"
  "s0700_f70_pL 0.70 0.01  0.012 0.0035 r0.02"
  "s0700_f95_pL 0.70 0.01  0.012 0.0004 r0.02"
  "s1000_f30_pL 1.00 0.01  0.012 0.0145 r0.02"
  "s1000_f95_pL 1.00 0.01  0.012 0.0003 r0.02"
  "s0300_f70_pH 0.30 0.01  0.012 0.0067 r0.1"
  "s0300_f95_pH 0.30 0.01  0.012 0.0008 r0.1"
  "s0700_f70_pH 0.70 0.01  0.012 0.0035 r0.1"
  "s0700_f95_pH 0.70 0.01  0.012 0.0004 r0.1"
)
DOSES=(${DOSES:-5 10 20})
set -u
i=$1; NSET=${2:-5}; NCH=${3:-3}
read -r lab irr ams f g pk <<< "${CASES[$i]}"
tag=as_f${f}_g${g}_a1
source stanage_handoff/v22_options.env
export XRSS=,as3up KEEPPH=1 XRS=,fhage
EXE=${AS_EXE:-scepter_paddy_fast_p}; K=${AS_K:-paddy/krxn_calibrated_v25_cal.in}
PROD=$(echo $PROD | sed "s#--k=[^ ]*#--k=$K#")
S="--fe3=$f --gt=$g --irralk=1 --irrsi=$irr --amsi=$ams --askd=paddy/arsenic/kd_constants_P$pk.json --as=${AS_SOIL:-10}:${AS_IRR:-100}:0.1 --exe=$EXE ${SOILX:-}"
OUT=${SCEPTER_OUT%/}/fac_$lab/; mkdir -p $OUT/paddy
export SCEPTER_OUT=$OUT
ST=${ST:-status/as_factorial_$lab.txt}
grep -q "^${tag}_spin rc 0" $ST 2>/dev/null || {
  python3 paddy/run_flooded_column.py ${tag}_spin fe2o2,amo2o,as3sp,as5sp,as5rd,as3ox,as5ra,as3oa 3 --fallow --basalt=0 $SOIL $P21 $S \
    > /dev/null 2>&1; echo "${tag}_spin rc $?" >> $ST; }
paddy/multiyear/run_chain.sh ${tag}_set $NSET 0 ${tag}_spin $ST $PROD $P21 $S --kaq=0.74 || exit 1
for d in 0 ${DOSES[@]}; do
  [ -d $OUT/paddy/flood_${tag}_st$d ] || cp -a $OUT/paddy/flood_${tag}_set_y${NSET}f $OUT/paddy/flood_${tag}_st$d
  NOPRUNE=1 paddy/multiyear/run_chain.sh ${tag}_d$d $NCH $d ${tag}_st$d $ST $PROD $P21 $S --kaq=0.74 &
done; wait
