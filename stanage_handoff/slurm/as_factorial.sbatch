#!/bin/bash
# Representative baseline Si x dose x f x phosphate subset (paddy/arsenic/run_as_factorial.sh): 14 soils, control + 5/10/20 t/ha.
# Build its executable under its own name first, so running jobs keep their binary:
#   module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
#   make --file=makefile_paddy_fast PROGRAM=scepter_paddy_fast_p OBJS=scepter_paddy_fast_p.o
#   mkdir -p logs status && sbatch --array=0-13 stanage_handoff/slurm/as_factorial.sbatch
#SBATCH --job-name=as_fact
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=16:00:00
#SBATCH --partition=sheffield
#SBATCH --account=default
#SBATCH --output=logs/as_factorial_%A_%a.out
module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # one BLAS thread per SCEPTER process: multithreaded OpenBLAS stalls basalt seasons (PSD step collapses dt)
python3 -c 'import numpy' || exit 1
[ -x scepter_paddy_fast_p ] || { echo "build scepter_paddy_fast_p first" >&2; exit 1; }
export SCEPTER_OUT=${SCEPTER_OUT:-/mnt/parscratch/users/$USER/scepter_as/}
bash paddy/arsenic/run_as_factorial.sh $SLURM_ARRAY_TASK_ID
