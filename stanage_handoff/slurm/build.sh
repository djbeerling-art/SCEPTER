#!/bin/bash
# Build the frozen v22 executable on Stanage (run once on a login node or in a short job).
# Module names are a guess: check with `module spider GCC` and `module spider OpenBLAS` and edit.
module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 2>/dev/null || module load GCC OpenBLAS
cd $HOME/SCEPTER && git fetch origin && git checkout claude/project-thread-lst1mi   # = pip6bj-pan + global-run prep
# v22 = scepter.f90 as of commit 162913a (the branch may later carry the v23 carbon fix)
git show 162913a:scepter.f90 > scepter.f90.v22 && cp scepter.f90 scepter.f90.head && cp scepter.f90.v22 scepter.f90
make --file=makefile_paddy PROGRAM=scepter_paddy_v22 OBJS=scepter_paddy_v22.o
mv scepter.f90.head scepter.f90
python3 -m pip install --user numpy scipy pandas scikit-learn SALib
