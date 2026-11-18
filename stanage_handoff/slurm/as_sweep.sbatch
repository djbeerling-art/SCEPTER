#!/bin/bash
# Arsenic response sweep: one array task per soil (27). From the repo root on the arsenic branch:
#   git fetch origin && git checkout claude/arsenic-sorption && git pull
#   module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
#   make --file=makefile_paddy_fast clean; make --file=makefile_paddy_fast
#   python3 -m pip install --user numpy   # once, if numpy is missing
#   mkdir -p logs status && sbatch --array=0-26 stanage_handoff/slurm/as_sweep.sbatch
# When all tasks finish:
#   python3 paddy/arsenic/as_sweep_summary.py $SCEPTER_OUT/paddy > as_sweep_summary.csv
# 5 settle + 3 control/basalt years per soil, ~3-5 core-hours; resubmitting resumes from the status file.
#SBATCH --job-name=as_sweep
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --partition=sheffield
#SBATCH --account=default
#SBATCH --output=logs/as_sweep_%A_%a.out
module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0   # same toolchain as the build
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # one BLAS thread per SCEPTER process: multithreaded OpenBLAS stalls basalt seasons (PSD step collapses dt)
export SCEPTER_OUT=${SCEPTER_OUT:-/mnt/parscratch/users/$USER/scepter_as/}
mkdir -p $SCEPTER_OUT/paddy logs status
python3 -c 'import numpy' || { echo 'numpy missing: python3 -m pip install --user numpy (with the modules above loaded)' >&2; exit 1; }
ST=status/as_sweep_$SLURM_ARRAY_TASK_ID.txt bash paddy/arsenic/run_as_sweep.sh $SLURM_ARRAY_TASK_ID ${NSET:-5} ${NCH:-3} ${DOSE:-20}
