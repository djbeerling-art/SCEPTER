#!/bin/bash
# Run the Morris screen job lists written by morris.py design: spin-ups, then controls and basalt runs.
# Usage (repository root, SCEPTER_OUT set): paddy/screen/run_screen.sh <dir> [parallel jobs, default 4]
# Rerunning skips jobs whose status line already reads rc 0.
set -u
export D=$(cd "$1" && pwd); NJ=${2:-4}
export RG=fe2o2,g2fe3,g2gr,g2so4,g2ch4,ch4o2,fe2gr,ch4pl,ch4eb,fe2fs,fesox,amo2o,g2n21,g2n2p,n2opl,g2gt,rtco2,g2o2,nh4up,rtoxy
export SOIL="--ch4 --n --kaol=0.05 --refsoil --gt=0.02 --nz=15 --exe=scepter_paddy_fast"
export SEAS="--season --nfert=120 --ch4ph=6.2:7.9:6.5 --nucl=cc:0.001 --kscale=maf:0.45,fsp:0.25"
job(){ tag=${1%% *}; o=${1#* }; grep -q "^$tag rc 0" "$D/status" 2>/dev/null && return
  case $tag in
    ms_*) python3 paddy/run_flooded_column.py $tag fe2o2,amo2o 3 --fallow --basalt=0 $SOIL $o > "$D/$tag.log" 2>&1;;
    mc_*) python3 paddy/run_flooded_column.py $tag $RG 0.31 $SOIL $SEAS $o > "$D/$tag.log" 2>&1;;
    *)    python3 paddy/run_flooded_column.py $tag $RG 0.31 $SOIL $SEAS --psd=data/psdrain_renforth_p80-100um.in $o > "$D/$tag.log" 2>&1;;
  esac
  rc=$?; r=$SCEPTER_OUT/paddy/flood_$tag
  [ "${tag:0:3}" != ms_ ] && find $r/prof -name "*.txt" -delete 2>/dev/null; find $r/flx -name "flx_psd*" -delete 2>/dev/null
  echo "$tag rc $rc" >> "$D/status"; }
export -f job
xargs -P $NJ -d '\n' -I{} bash -c 'job "{}"' < "$D/jobs_spin.txt"
cat "$D/jobs_ctrl.txt" "$D/jobs_dose.txt" | xargs -P $NJ -d '\n' -I{} bash -c 'job "{}"'
echo done >> "$D/status"
