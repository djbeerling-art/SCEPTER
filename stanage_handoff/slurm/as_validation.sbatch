#!/bin/bash
# Bengal validation clusters (paddy/arsenic/run_as_validation.sh): 6 clusters x 3 f levels, control + 3 Rajmahal doses.
# Build its executable under its own name first, so running jobs keep their binary:
#   module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
#   make --file=makefile_paddy_fast PROGRAM=scepter_paddy_fast_p OBJS=scepter_paddy_fast_p.o
#   mkdir -p logs status && sbatch --array=0-17 stanage_handoff/slurm/as_validation.sbatch
#SBATCH --job-name=as_val
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=16:00:00
#SBATCH --partition=sheffield
#SBATCH --account=default
#SBATCH --output=logs/as_validation_%A_%a.out
module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # one BLAS thread per SCEPTER process: multithreaded OpenBLAS stalls basalt seasons (PSD step collapses dt)
python3 -c 'import numpy' || exit 1
[ -x scepter_paddy_fast_p ] || { echo "build scepter_paddy_fast_p first" >&2; exit 1; }
export SCEPTER_OUT=${SCEPTER_OUT:-/mnt/parscratch/users/$USER/scepter_as/}
bash paddy/arsenic/run_as_validation.sh $SLURM_ARRAY_TASK_ID
