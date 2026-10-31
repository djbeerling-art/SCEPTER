#!/bin/bash
# One Slurm array task per archetype. Submit from the repo root:
#   mkdir -p logs status    # once (build.sh does it): Slurm will not create logs/ and the task dies without a log
#   sbatch --array=1-201%100 stanage_handoff/slurm/arch_array.sbatch stanage_handoff/arch/archetypes_v22.csv
# DOSES, YEARS, STATUS_DIR (default status) and SCEPTER_OUT can be overridden for test runs kept apart from production,
# but ONLY via --export=ALL,...: Stanage sets SBATCH_EXPORT=HOME,PATH,TERM, so inline VAR=x sbatch ... is silently dropped
# and the task runs with production defaults. Example (1-yr build check of irr061 = row 4 of pilot10.csv):
#   sbatch --export=ALL,DOSES="0 20",YEARS=1,STATUS_DIR=status_check,SCEPTER_OUT=/mnt/parscratch/users/$USER/scepter_check/ \
#          --array=4 --time=04:00:00 stanage_handoff/slurm/arch_array.sbatch stanage_handoff/arch/pilot10.csv
# 10-yr chains at 0/10/20 t/ha: ~8-9 core-hours per archetype in the cloud container (basalt years slow to ~25 min
# after year 5), so 24 h wall time leaves margin; the script resumes if resubmitted.
#SBATCH --job-name=rice23
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --partition=sheffield
#SBATCH --account=default
#SBATCH --output=logs/rice23_%A_%a.out
source stanage_handoff/slurm/modules.sh
source stanage_handoff/v22_options.env
export SCEPTER_OUT=${SCEPTER_OUT:-/mnt/parscratch/users/$USER/scepter/}
STATUS_DIR=${STATUS_DIR:-status}
mkdir -p $SCEPTER_OUT/paddy $STATUS_DIR
row=$(python3 -c "
import csv,sys; r=list(csv.DictReader(open('$1')))[$SLURM_ARRAY_TASK_ID-1]
a=r['args']+(' --awd' if r['regime']=='rf' else ''); print(r['arch']+'|'+a)")
IFS='|' read a args <<< "$row"
bash stanage_handoff/run_arch_v22.sh $a "$args" "${DOSES:-0 10 20}" ${YEARS:-10} $STATUS_DIR/$a.txt
