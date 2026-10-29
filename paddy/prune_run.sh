#!/bin/bash
# Shrink finished run directories to what paddy/ghg_budget.py reads (season runs) or what a restart reads (spin-ups).
# Usage: paddy/prune_run.sh season|spin <run dir> ...
# season: keep flx_gas-*, flx_co2sp-ALK, flx_sld-*, flx_aq-{ca,mg,na,k,fe2}, zint_* (layer budgets), int_flx_ Ca/Mg/Fe2/ALK/solids (exact season totals); drop prof/ and everything else in flx/.
# spin:   drop flx/ entirely; prof/ (restart state) is kept.
# KEEPPH=1 (environment) keeps the season's porewater profiles prof_aq-NNN.txt (pH time series); KEEPALK=1 also keeps
# prof_aq(tot)-, prof_aq(ads)-, prof_gas(aq)- and prof_sld-NNN.txt (paddy/alk_flux.py).
mode=$1; shift
for r in "$@"; do
  if [ "$mode" = spin ]; then find "$r/flx" -type f -delete 2>/dev/null; continue; fi
  if [ -n "${KEEPALK:-}" ]; then find "$r/prof" -type f ! -name 'prof_aq-[0-9]*.txt' ! -name 'prof_aq(tot)-[0-9]*.txt' ! -name 'prof_gas(aq)-[0-9]*.txt' ! -name 'prof_sld-[0-9]*.txt' ! -name 'prof_aq(ads)-[0-9]*.txt' -delete 2>/dev/null
  elif [ -n "${KEEPPH:-}" ]; then find "$r/prof" -type f ! -name 'prof_aq-[0-9]*.txt' -delete 2>/dev/null
  else find "$r/prof" -type f -delete 2>/dev/null; fi
  find "$r/flx" -type f ! -name 'flx_gas-*' ! -name 'flx_co2sp-ALK.txt' ! -name 'flx_sld-*' \
    ! -name 'flx_aq-ca.txt' ! -name 'flx_aq-mg.txt' ! -name 'flx_aq-na.txt' ! -name 'flx_aq-k.txt' ! -name 'flx_aq-fe2.txt' ! -name 'zint_*' ! -name 'int_flx_aq-ca.txt' ! -name 'int_flx_aq-mg.txt' ! -name 'int_flx_aq-fe2.txt' ! -name 'int_flx_sld-*' ! -name 'int_flx_co2sp-ALK.txt' -delete 2>/dev/null
done
