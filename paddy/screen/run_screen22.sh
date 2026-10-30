#!/bin/bash
# Run the v22 Morris screen job lists written by morris22.py design: spin-ups, 7-yr settles, then one-year controls and
# basalt runs (season + fallow, run_chain.sh). Frozen v22 production options as the archetype runs (runZa.sh).
# Usage (repository root, SCEPTER_OUT set): paddy/screen/run_screen22.sh <dir> [parallel jobs, default 4]. Resumable.
set -u
export D=$(cd "$1" && pwd); NJ=${2:-4}; export ST=$D/status
source <(sed -n '1,/^export -f spinset season/p' /home/user/runs/pan/runV.sh | sed "s#ST=/home/user/runs/pan/statusD#ST=$ST#")
export ST=$D/status; unset KEEPALK
export PROD="$(echo $PROD | sed 's/scepter_paddy_v21/scepter_paddy_v22/') --rtsi=2e-3:2.7e-4:0.2"
export SOIL="$(echo $SOIL | sed 's/scepter_paddy_v2[01]/scepter_paddy_v22/')"
job(){ set -- $1; tag=$1; shift
  case $tag in
    ms22_*) grep -q "^$tag rc 0" $ST && return
      python3 paddy/run_flooded_column.py $tag fe2o2,amo2o 3 --fallow --basalt=0 $SOIL $P21 "$@" > $D/$tag.log 2>&1
      echo "$tag rc $?" >> $ST; paddy/prune_run.sh spin $SCEPTER_OUT/paddy/flood_$tag;;
    mv22_*) sp=$1; shift; paddy/multiyear/run_chain.sh $tag 7 0 $sp $ST $PROD $P21 "$@" &&
      for y in 1 2 3 4 5 6; do rm -rf $SCEPTER_OUT/paddy/flood_${tag}_y${y}s $SCEPTER_OUT/paddy/flood_${tag}_y${y}f; done;;
    *) grep -q "^${tag}_y1f rc 0" $ST && return; st=$1; d=$2; shift 2; o=$SCEPTER_OUT/paddy/flood_r$tag; [ -d $o ] || cp -a $SCEPTER_OUT/paddy/flood_${st}_y7f $o
      paddy/multiyear/run_chain.sh $tag 1 $d r$tag $ST $PROD $P21 "$@" && { rm -rf $o
        paddy/prune_run.sh season $SCEPTER_OUT/paddy/flood_${tag}_y1s; paddy/prune_run.sh season $SCEPTER_OUT/paddy/flood_${tag}_y1f; };;
  esac; }
export -f job
for p in spin setl; do xargs -P $NJ -d '\n' -I{} bash -c 'job "{}"' < $D/jobs_$p.txt; done
# spin-ups are no longer needed once every settle has run (disk)
grep -q '^setl done' $ST || { [ $(grep -c '_y7f rc 0' $ST) -ge $(wc -l < $D/jobs_setl.txt) ] && { for t in $(cut -d' ' -f1 $D/jobs_spin.txt); do rm -rf $SCEPTER_OUT/paddy/flood_$t; done; echo 'setl done' >> $ST; }; }
cat $D/jobs_ctrl.txt $D/jobs_dose.txt | xargs -P $NJ -d '\n' -I{} bash -c 'job "{}"'
echo "screen done" >> $ST
