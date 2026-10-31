#!/bin/bash
# One archetype, frozen v22 protocol (identical to the GMD archetype runs, runR.sh):
#   3-yr drained spin-up (gs22_<a>) -> 7-yr settle without rock (st22_<a>, kept as y7f) -> per dose, a copy of the settled
#   state and an N-yr chain (g22_<a>_<dose>) of 113-d season + 252-d fallow, basalt applied at the start of every season.
# Usage (repository root, SCEPTER_OUT set, v22_options.env sourced):
#   run_arch_v22.sh <arch> "<driver args>" "<doses e.g. 0 10 20>" <years> <status file>
# <driver args> come from archetypes.csv 'args' (rock, irralk, soilcc, cec, temp, om, fe3) plus --awd for rainfed archetypes.
# Rerunning resumes: finished phases are skipped via the status file. Disk: ~30 MB per pruned year.
set -u
a=$1; args=$2; doses=$3; ny=$4; ST=$5
ok(){ grep -q "^$1 rc 0" "$ST" 2>/dev/null; }
# 3-yr drained spin-up (--fallow; args passed exactly as in the GMD runs)
if ! ok gs22_$a; then
  python3 paddy/run_flooded_column.py gs22_$a fe2o2,amo2o 3 --fallow --basalt=0 $SOIL $P21 $args > /dev/null 2>&1
  echo "gs22_$a rc $?" >> "$ST"; ok gs22_$a || exit 1
  paddy/prune_run.sh spin $SCEPTER_OUT/paddy/flood_gs22_$a
fi
paddy/multiyear/run_chain.sh st22_$a 7 0 gs22_$a "$ST" $PROD $P21 $args || exit 1
for y in 1 2 3 4 5 6; do rm -rf $SCEPTER_OUT/paddy/flood_st22_${a}_y${y}[sf]; done
for d in $doses; do
  o=$SCEPTER_OUT/paddy/flood_gx22_${a}_$d
  ok g22_${a}_${d}_y${ny}f && continue
  [ -d $o ] || cp -a $SCEPTER_OUT/paddy/flood_st22_${a}_y7f $o
  paddy/multiyear/run_chain.sh g22_${a}_$d $ny $d gx22_${a}_$d "$ST" $PROD $P21 $args || exit 1
  rm -rf $o
done
