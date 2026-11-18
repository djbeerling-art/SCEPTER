#!/bin/bash
# Bengal validation clusters (David, 3 Nov 2026): full SCEPTER at six soil clusters chosen from the Chowdhury et al. 2017
# paddy-soil survey (1,209 georeferenced soils), spanning irrigation source, soil As, pH and geology. Used to test the
# correction-A Bengal map (bengal_v1) out of sample. Each task = one cluster x one ferrihydrite-share level, taken through
# the sweep workflow (drained spin-up, 5-yr settle without rock), then a control and three Rajmahal basalt doses from the
# same settled state (year-matched control).
#
# Usage (repository root, SCEPTER_OUT set): paddy/arsenic/run_as_validation.sh <task 0-17> [settle yrs 5] [chain yrs 3]
#   task = 3*cluster + f level (0 low ~0.5, 1 central ~0.8, 2 high ~0.9; goethite ladder as realised in the factorial)
# Output: $SCEPTER_OUT/val_<cluster>_f<level>/paddy; summarise with as_sweep_summary.py <dir> <dose>.
#
# Inputs per cluster (upazila medians of the survey; private data, values only):
#   ferrihydrite = 0.012 x total Fe / 26,800 mg/kg (survey median), goethite ladder scaled the same way (keeps f);
#   labile soil As = 0.15 x total As (oxalate/total As ~0.13-0.21 in Stroud et al. 2011);
#   irrigation Si: groundwater clusters 1.8 mM (in the factorial, settled porewater Si sat far below irrigation Si:
#     1.0 mM gave ~430 uM; 1.8 mM aims at the ~600 uM measured by Roberts et al. 2011), surface-water clusters 0.15 mM
#     (Mekong mean SiO2 8.85 mg/L; Ganges-Brahmaputra floodwater Si unmeasured, assumed similar);
#   irrigation As: groundwater clusters provisional local values (Sreenagar well 400 ug/L, Roberts et al. 2011),
#     surface water 2 ug/L; the uptake ratio is insensitive to this (As is trace);
#   pH via irrigation alkalinity: irralk 0 (measured pH ~4.7-5.1), 1 (5.6-6.4), 5 + 2% calcite (Satkhira, 7.3);
#   phosphate low (P/Fe 0.02 refit; high P strengthens suppression, so low P is the cautious case);
#   doses: groundwater clusters 5/10/20 t/ha; surface-water (low-Si) clusters 2/5/10 t/ha (low-dose risk, boundary test).
# Cluster: label irrSi(mM) fe3 As_soil(mg/kg) As_irr(ug/L) irralk soilcc doses
CL=(
  "munshiganj_gw 1.8 0.0142 1.18 400 1 0    5,10,20"
  "comilla_gw    1.8 0.0121 1.39 200 1 0    5,10,20"
  "satkhira_gw   1.8 0.0158 2.15 150 5 0.02 5,10,20"
  "madhupur_gw   1.8 0.0059 0.39 10  0 0    5,10,20"
  "gopalganj_sw  0.15 0.0139 1.48 2  1 0    2,5,10"
  "lakshmipur_sw 0.15 0.0112 0.47 2  0 0    2,5,10"
)
GTF=(2.5 0.58 0.083)   # goethite / ferrihydrite ratio for f ~0.5 / 0.8 / 0.9 (factorial: gt 0.03/0.007/0.001 at fe3 0.012)
set -u
i=$1; NSET=${2:-5}; NCH=${3:-3}
read -r lab irr f3 aso air alk scc dl <<< "${CL[$((i/3))]}"
[ -n "${IRR:-}" ] && { irr=$IRR; lab=${lab}_si$IRR; }
[ -n "${ASIRR:-}" ] && { air=$ASIRR; lab=${lab}_as$ASIRR; }   # ASIRR (env, ug/L) overrides irrigation As (irrigation-As sweep)   # IRR (env, mM) overrides irrigation Si; separate output/status
k=$((i%3)); g=$(python3 -c "print(round($f3*${GTF[$k]},5))")
DOSES=(${DOSES:-${dl//,/ }})
tag=as_f${f3}_g${g}_a${alk}
source stanage_handoff/v22_options.env
export XRSS=,as3up KEEPPH=1 XRS=,fhage
EXE=${AS_EXE:-scepter_paddy_fast_p}; K=${AS_K:-paddy/krxn_calibrated_v25_cal.in}
PROD=$(echo $PROD | sed "s#--k=[^ ]*#--k=$K#")
S="--fe3=$f3 --gt=$g --irralk=$alk --irrsi=$irr --amsi=0.01 --soilcc=$scc --askd=paddy/arsenic/kd_constants_Pr0.02.json --as=$aso:$air:0.1 --rock=rajmahal --exe=$EXE ${SOILX:-}"
OUT=${SCEPTER_OUT%/}/val_${lab}_f$k/; mkdir -p $OUT/paddy
export SCEPTER_OUT=$OUT
ST=${ST:-status/as_validation_${lab}_f$k.txt}
grep -q "^${tag}_spin rc 0" $ST 2>/dev/null || {
  python3 paddy/run_flooded_column.py ${tag}_spin fe2o2,amo2o,as3sp,as5sp,as5rd,as3ox,as5ra,as3oa 3 --fallow --basalt=0 $SOIL $P21 $S \
    > /dev/null 2>&1; echo "${tag}_spin rc $?" >> $ST; }
paddy/multiyear/run_chain.sh ${tag}_set $NSET 0 ${tag}_spin $ST $PROD $P21 $S --kaq=0.74 || exit 1
for d in 0 ${DOSES[@]}; do
  [ -d $OUT/paddy/flood_${tag}_st$d ] || cp -a $OUT/paddy/flood_${tag}_set_y${NSET}f $OUT/paddy/flood_${tag}_st$d
  NOPRUNE=1 paddy/multiyear/run_chain.sh ${tag}_d$d $NCH $d ${tag}_st$d $ST $PROD $P21 $S --kaq=0.74 &
done; wait
