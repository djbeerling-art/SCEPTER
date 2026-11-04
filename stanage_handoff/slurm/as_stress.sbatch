#!/bin/bash
# Arsenic stress test: 4 soils nearest the suppress/mobilise boundary x 5 scenarios = 20 tasks, one year each.
#   module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
#   mkdir -p logs status && sbatch --array=0-19 stanage_handoff/slurm/as_stress.sbatch
# Then per scenario: python3 paddy/arsenic/as_sweep_summary.py $SCEPTER_OUT/stress_<sc>/paddy <dose>
#SBATCH --job-name=as_stress
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --partition=sheffield
#SBATCH --account=default
#SBATCH --output=logs/as_stress_%A_%a.out
module purge; module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # one BLAS thread per SCEPTER process: multithreaded OpenBLAS stalls basalt seasons (PSD step collapses dt)
python3 -c 'import numpy' || exit 1
export SCEPTER_OUT=${SCEPTER_OUT:-/mnt/parscratch/users/$USER/scepter_as/}
SOILS=(10 11 19 20); SC=(${SCN:-dose5 dose10 dose40 ox01 fsite2}); n=${#SC[@]}   # SCN=worst40: sbatch --export=ALL,SCN=worst40 --array=0-3
t=$SLURM_ARRAY_TASK_ID
bash paddy/arsenic/run_as_stress.sh ${SOILS[$((t/n))]} ${SC[$((t%n))]}
