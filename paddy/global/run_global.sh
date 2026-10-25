#!/bin/bash
# Run every archetype of archetypes.csv (archetypes.py): a 3-yr fallow spin-up, then one season at each basalt dose.
# Usage: paddy/global/run_global.sh <outdir with archetypes.csv> [doses, default "0 10 20"] [parallel jobs, default 4]
# Run from the repository root with SCEPTER_OUT set; logs and status go to <outdir>. Tags: gs_<arch> (spin-up),
# g_<arch>_<dose> (season). Rerunning skips jobs whose status line already reads rc 0.
set -u
export D=$(cd "$1" && pwd); export DOSES=${2:-"0 10 20"}; NJ=${3:-4}
export RG=fe2o2,g2fe3,g2gr,g2so4,g2ch4,ch4o2,fe2gr,ch4pl,ch4eb,fe2fs,fesox,amo2o,g2n21,g2n2p,n2opl,g2gt,rtco2,g2o2,nh4up,rtoxy
export SOIL="--ch4 --n --fe3=0.012 --kaol=0.05 --refsoil --gt=0.02 --nz=15 --exe=scepter_paddy_fast"
export SEAS="--k=paddy/krxn_calibrated_v5.in --season --nfert=120 --ch4ph=6.2:7.9:6.5 --nucl=cc:0.001 --kscale=maf:0.45,fsp:0.25"
export PSD=data/psdrain_renforth_p80-100um.in
done_ok(){ grep -q "^$1 rc 0" "$D/status" 2>/dev/null; }
arch_job(){ IFS=, read -r aid regime args <<< "$1"
  [ "$regime" = irr ] && w=--phd=14 || w=--awd
  if ! done_ok gs_$aid; then
    python3 paddy/run_flooded_column.py gs_$aid fe2o2,amo2o 3 --fallow --basalt=0 $SOIL $args > "$D/gs_$aid.log" 2>&1
    echo "gs_$aid rc $?" >> "$D/status"; done_ok gs_$aid || return
  fi
  for d in $DOSES; do tag=g_${aid}_$d; done_ok $tag && continue
    p=""; [ "$d" != 0 ] && p=--psd=$PSD
    python3 paddy/run_flooded_column.py $tag $RG 0.31 --rstrt=gs_$aid --basalt=$d $SOIL $SEAS $w $p $args > "$D/$tag.log" 2>&1
    rc=$?; [ $rc = 0 ] && [ "${tag:0:3}" != ms_ ] && ! awk 'END{exit !($1+0>=0.3095)}' $SCEPTER_OUT/paddy/flood_$tag/flx/flx_gas-pch4.txt 2>/dev/null && rc=98   # stopped before season end
    o=$SCEPTER_OUT/paddy/flood_$tag
    find $o/prof -name "*.txt" -delete 2>/dev/null; find $o/flx -name "flx_psd*" -delete 2>/dev/null   # keep disk use ~3 MB per season run
    echo "$tag rc $rc" >> "$D/status"
  done; }
export -f arch_job done_ok
python3 -c "
import csv,sys
for r in csv.DictReader(open('$D/archetypes.csv')): print(r['arch']+','+r['regime']+','+r['args'])" \
  | xargs -P $NJ -d '\n' -I{} bash -c 'arch_job "{}"'
echo done >> "$D/status"
