#!/bin/bash
# Build the production executable on Stanage (run once on a login node or in a short interactive job).
# v23 (David, 31 Oct 2026: v23 replaces v22) = scepter.f90 at commit 267a856: the v22 model plus the carbon-closure fix
# (oxalate/acetate C from g2 decay returned to CO2 when oxa/ac are untracked). Same options and krxn_calibrated_v22.in.
# Module names: check with `module spider GCC` / `module spider OpenBLAS` and edit if these fail.
module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 2>/dev/null || module load GCC OpenBLAS
cd $HOME/SCEPTER && git fetch origin && git checkout claude/project-thread-lst1mi && git pull --ff-only
git show 267a856:scepter.f90 > scepter.f90.v23 && cmp -s scepter.f90 scepter.f90.v23 || echo "WARNING: scepter.f90 differs from v23 (267a856)"
cp scepter.f90 scepter.f90.head && cp scepter.f90.v23 scepter.f90
make --file=makefile_paddy PROGRAM=scepter_paddy_v23 OBJS=scepter_paddy_v23.o
mv scepter.f90.head scepter.f90
python3 -m pip install --user numpy scipy pandas scikit-learn SALib matplotlib
