#!/bin/bash
# Build the production executable on Stanage (run once on a login node or in a short interactive job).
# v23 (David, 31 Oct 2026: v23 replaces v22) = scepter.f90 at commit 267a856: the v22 model plus the carbon-closure fix
# (oxalate/acetate C from g2 decay returned to CO2 when oxa/ac are untracked). Same options and krxn_calibrated_v22.in.
# Modules come from modules.sh (shared with arch_array.sbatch).
set -e
[ -d $HOME/SCEPTER ] || git clone https://github.com/djbeerling-art/SCEPTER $HOME/SCEPTER
cd $HOME/SCEPTER && git fetch origin && git checkout claude/project-thread-lst1mi && git pull --ff-only   # = pip6bj-pan + global-run prep
source stanage_handoff/slurm/modules.sh
git show 267a856:scepter.f90 > scepter.f90.v23
cmp -s scepter.f90 scepter.f90.v23 || echo "NOTE: branch scepter.f90 differs from v23 (267a856); building v23"
cp scepter.f90 scepter.f90.head && cp scepter.f90.v23 scepter.f90
# restore the branch's scepter.f90 even if make fails, so the next git checkout/pull is not blocked
make --file=makefile_paddy PROGRAM=scepter_paddy_v23 OBJS=scepter_paddy_v23.o || { mv scepter.f90.head scepter.f90; exit 1; }
mv scepter.f90.head scepter.f90
python3 -m pip install --user numpy scipy pandas scikit-learn SALib matplotlib
# Slurm does not create the --output directory: a job whose logs/ is missing fails at once with no log at all.
mkdir -p logs status
ls -l scepter_paddy_v23 && module list
