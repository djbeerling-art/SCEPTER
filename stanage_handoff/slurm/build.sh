#!/bin/bash
# Build the frozen v22 executable on Stanage (run once on a login node or in a short job).
# Modules come from modules.sh (shared with arch_array.sbatch).
set -e
[ -d $HOME/SCEPTER ] || git clone https://github.com/djbeerling-art/SCEPTER $HOME/SCEPTER
cd $HOME/SCEPTER && git fetch origin && git checkout claude/project-thread-lst1mi && git pull --ff-only   # = pip6bj-pan + global-run prep
source stanage_handoff/slurm/modules.sh
# v22 = scepter.f90 as of commit 162913a (the branch may later carry the v23 carbon fix)
git show 162913a:scepter.f90 > scepter.f90.v22 && cp scepter.f90 scepter.f90.head && cp scepter.f90.v22 scepter.f90
# restore the branch's scepter.f90 even if make fails, so the next git checkout/pull is not blocked
make --file=makefile_paddy PROGRAM=scepter_paddy_v22 OBJS=scepter_paddy_v22.o || { mv scepter.f90.head scepter.f90; exit 1; }
mv scepter.f90.head scepter.f90
python3 -m pip install --user numpy scipy pandas scikit-learn SALib
# Slurm does not create the --output directory: a job whose logs/ is missing fails at once with no log at all.
mkdir -p logs status
ls -l scepter_paddy_v22 && module list
