#!/bin/bash
# One module set for the build AND the array jobs, so the executable runs against the libgfortran/OpenBLAS
# it was linked with, and python3 is the same interpreter the packages were pip-installed into.
# EasyBuild 2023a names (Stanage). Unverified: if a load fails, run `module spider GCC OpenBLAS Python` and edit.
module purge 2>/dev/null
module load GCC/12.3.0 OpenBLAS/0.3.23-GCC-12.3.0 Python/3.11.3-GCCcore-12.3.0 \
  || { echo "modules.sh: versioned load failed, falling back to defaults" >&2; module load GCC OpenBLAS Python; }
