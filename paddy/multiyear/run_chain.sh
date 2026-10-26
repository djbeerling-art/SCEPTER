#!/bin/bash
# Multi-year paddy run as a chain of restarts: each year = one 113-d rice season (basalt applied at its start if dose > 0,
# urea, plants) followed by a drained, crop-free fallow (water table 0.5 m, no plant reactions, no fertiliser).
# Usage (repository root, SCEPTER_OUT set): paddy/multiyear/run_chain.sh <name> <years> <dose t/ha> <spin-up tag> <status file> [driver options]
# (no methanogen pH response unless --ch4ph=... is passed)
# Tags: <name>_y<n>s (season), <name>_y<n>f (fallow). Rerunning skips years already finished (status line rc 0).
# The annual OM input (--om=, g C m-2 yr-1, default 450) is split between the phases: a fraction OMF (environment,
# default 0.7: roots, exudates and residue incorporated at land preparation) during the season, the rest (stubble,
# weeds) during the fallow; FOPT (environment) adds driver options to the fallow only (e.g. "--temp=18 --fsat=0.3"); each phase gets the matching rate (g C m-2 yr-1) for its length.
set -u
name=$1; ny=$2; dose=$3; prev=$4; st=$5; shift 5
om=$(printf '%s\n' "$@" | grep '^--om=' | cut -d= -f2); om=${om:-450}; OMF=${OMF:-0.7}
oms=$(python3 -c "print(round($om*$OMF/0.31,2))"); omf=$(python3 -c "print(round($om*(1-$OMF)/0.69,2))")
set -- $(printf '%s\n' "$@" | grep -v '^--om=')
RS=fe2o2,g2fe3,g2gr,g2so4,g2ch4,ch4o2,fe2gr,ch4pl,ch4eb,fe2fs,fesox,amo2o,g2n21,g2n2p,n2opl,g2gt,rtco2,g2o2,nh4up,rtoxy,groxy
RF=fe2o2,g2fe3,g2gr,g2so4,g2ch4,ch4o2,fe2gr,ch4eb,fe2fs,fesox,amo2o,g2n21,g2n2p,g2gt,g2o2,groxy
C="--ch4 --n --kaol=0.05 --refsoil --gt=0.02 --nz=15 --exe=scepter_paddy_fast --nucl=cc:0.001 --kscale=maf:0.45,fsp:0.25"
ok(){ grep -q "^$1 rc 0" "$st" 2>/dev/null; }
for y in $(seq 1 $ny); do
  ts=${name}_y${y}s; tf=${name}_y${y}f
  if ! ok $ts; then
    p=""; [ "$dose" != 0 ] && p=--psd=data/psdrain_renforth_p80-100um.in
    python3 paddy/run_flooded_column.py $ts $RS 0.31 --season --nfert=120 --rstrt=$prev --basalt=$dose --om=$oms $p $C "$@" > /dev/null 2>&1
    rc=$?; awk 'END{exit !($1+0>=0.3095)}' $SCEPTER_OUT/paddy/flood_$ts/flx/flx_gas-pch4.txt 2>/dev/null || rc=98
    find $SCEPTER_OUT/paddy/flood_$ts/flx -name 'flx_psd*' -delete 2>/dev/null
    echo "$ts rc $rc" >> "$st"; [ $rc = 0 ] || exit 1
  fi
  if ! ok $tf; then
    # the fallow takes the same soil options but no water-management flags (drained by --fallow)
    python3 paddy/run_flooded_column.py $tf $RF 0.69 --fallow --rstrt=$ts --basalt=0 --om=$omf $C $(printf '%s\n' "$@" | grep -v -e '^--awd$' -e '^--phd' -e '^--k=') --k=$(printf '%s\n' "$@" | grep '^--k=' | cut -d= -f2) ${FOPT:-} > /dev/null 2>&1
    rc=$?; awk 'END{exit !($1+0>=0.6895)}' $SCEPTER_OUT/paddy/flood_$tf/flx/flx_gas-pch4.txt 2>/dev/null || rc=98
    find $SCEPTER_OUT/paddy/flood_$tf/flx -name 'flx_psd*' -delete 2>/dev/null
    [ $y -gt 1 ] && find $SCEPTER_OUT/paddy/flood_${name}_y$((y-1))[sf]/prof -type f -delete 2>/dev/null
    echo "$tf rc $rc" >> "$st"; [ $rc = 0 ] || exit 1
  fi
  [ -d $SCEPTER_OUT/paddy/flood_$ts/prof ] && find $SCEPTER_OUT/paddy/flood_$ts/prof -type f -delete
  prev=$tf
done
