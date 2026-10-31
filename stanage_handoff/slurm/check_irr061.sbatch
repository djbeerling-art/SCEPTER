#!/bin/bash
# Build check: irr061, full protocol (3-yr spin-up, 7-yr settle) then 1 yr at 0 and 20 t/ha, compared with the cloud
# v23 reference. Submit from the repo root after build.sh:  sbatch stanage_handoff/slurm/check_irr061.sbatch
# Expected ~1.5-2 h on one core. Result: tail logs/check_irr061_<jobid>.out (ends PASS or lists failures).
#SBATCH --job-name=chk23
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --partition=sheffield
#SBATCH --account=default
#SBATCH --output=logs/check_irr061_%j.out
source stanage_handoff/slurm/modules.sh
source stanage_handoff/v22_options.env
export SCEPTER_OUT=${SCEPTER_OUT:-/mnt/parscratch/users/$USER/scepter_check/}
mkdir -p $SCEPTER_OUT/paddy status_check
t0=$(date +%s)
bash stanage_handoff/run_arch_v22.sh irr061 "--rock=glb04 --irralk=1.6 --cec=16.9 --temp=28.7 --om=506 --fe3=0.0112" "0 20" 1 status_check/irr061.txt
echo "rc $? after $(( ($(date +%s)-t0)/60 )) min"; cat status_check/irr061.txt
python3 stanage_handoff/compare_ref.py $SCEPTER_OUT/paddy
