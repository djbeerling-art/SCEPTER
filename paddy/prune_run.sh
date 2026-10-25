#!/bin/bash
# Shrink finished run directories to what paddy/ghg_budget.py reads (season runs) or what a restart reads (spin-ups).
# Usage: paddy/prune_run.sh season|spin <run dir> ...
# season: keep flx_gas-*, flx_co2sp-ALK, flx_sld-*, flx_aq-{ca,mg,na,k}; drop prof/ and everything else in flx/.
# spin:   drop flx/ entirely; prof/ (restart state) is kept.
mode=$1; shift
for r in "$@"; do
  if [ "$mode" = spin ]; then find "$r/flx" -type f -delete 2>/dev/null; continue; fi
  find "$r/prof" -type f -delete 2>/dev/null
  find "$r/flx" -type f ! -name 'flx_gas-*' ! -name 'flx_co2sp-ALK.txt' ! -name 'flx_sld-*' \
    ! -name 'flx_aq-ca.txt' ! -name 'flx_aq-mg.txt' ! -name 'flx_aq-na.txt' ! -name 'flx_aq-k.txt' -delete 2>/dev/null
done
