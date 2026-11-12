"""Add amorphous silica (silica gel) to a finished run's restart state, for the Si-only arsenic test from settled soils.
Usage: python3 paddy/arsenic/add_silica.py <run dir> <wt fraction of solids> [depth m, default 0.2]
Adds amsi (mol m-3) over z < depth and removes the same solid volume of inert (kaolinite molar volume), so porosity is
unchanged. Rewrites <run dir>/prof/prof_sld-save.txt; keeps the original as prof_sld-save.nosi.txt.
The run that restarts from it should use --kscale=...,amsi:30 (as v22/v25 production) so the gel dissolves in weeks."""
import sys, os, shutil
import numpy as np
MV_AMSI, MV_INRT, MW_AMSI = 25.739, 99.52, 60.085   # cm3/mol, cm3/mol (kaolinite), g/mol (scepter.f90)
d, wt = sys.argv[1], float(sys.argv[2]); depth = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2
f = os.path.join(d, 'prof', 'prof_sld-save.txt'); fb = os.path.join(d, 'prof', 'bsd-save.txt')
lines = open(f).read().splitlines(); head = lines[0]
names = ['z'] + [head[4+5*i:9+5*i].strip() for i in range((len(head.rstrip())-4)//5)]
x = np.array([[float(v) for v in l.split()] for l in lines[1:] if l.strip()])
hb = open(fb).readline().split(); b = np.loadtxt(fb, skiprows=1)
poro, dens = b[:, hb.index('poro')], b[:, hb.index('dens[g/cm3]')]
if 'amsi' not in names or 'inrt' not in names: sys.exit(f'add_silica.py: no amsi/inrt in {names}')
top = x[:, 0] < depth
add = wt*(1-poro)*dens*1e6/MW_AMSI                          # mol m-3 bulk
x[top, names.index('amsi')] += add[top]
x[top, names.index('inrt')] = np.maximum(x[top, names.index('inrt')] - add[top]*MV_AMSI/MV_INRT, 0)
shutil.copy(f, f.replace('.txt', '.nosi.txt'))
with open(f, 'w') as o:
    o.write(head+'\n')
    for r in x: o.write(''.join('%27.16E' % v for v in r)+'\n')
print(f'added {add[top].mean():.0f} mol/m3 amsi over z<{depth} m')
