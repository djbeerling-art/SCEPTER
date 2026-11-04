#!/bin/bash
# Si x f boundary test: 8 soils (baseline Si ~50/100/300/220 uM x two ferrihydrite shares), control + 1/4/20 t/ha basalt.
#   module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
#   mkdir -p logs status && sbatch --array=0-7 stanage_handoff/slurm/si_boundary.sbatch
# Then per case and dose: python3 paddy/arsenic/as_sweep_summary.py $SCEPTER_OUT/sib_<label>/paddy <dose>
#SBATCH --job-name=si_bound
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=16:00:00
#SBATCH --partition=sheffield
#SBATCH --account=default
#SBATCH --output=logs/si_boundary_%A_%a.out
module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # one BLAS thread per SCEPTER process: multithreaded OpenBLAS stalls basalt seasons (PSD step collapses dt)
python3 -c 'import numpy' || exit 1
export SCEPTER_OUT=${SCEPTER_OUT:-/mnt/parscratch/users/$USER/scepter_as/}
bash paddy/arsenic/run_si_boundary.sh $SLURM_ARRAY_TASK_ID
