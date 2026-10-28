# Puddling/tillage before a rice season: homogenise the plough layer of a finished run's solid restart profile.
# Usage: python3 paddy/multiyear/till.py <run dir> [depth m, default 0.16]
# Mixes the soil solids (OM g2, Fe(III) oxides amfe3/gt, green rust, siderite, FeS, calcite, clays, amorphous silica,
# urea) over cells with z < depth, conserving each one's inventory (uniform grid); rock minerals are left alone because
# their particle-size distributions are tracked per depth. Without it, Fe(II) diffusing up and oxidising at the surface
# piles Fe(III) into the top cell and strips it from the rooted layer below over successive seasons.
# Rewrites <run dir>/prof/prof_sld-save.txt (the file a restart reads) and keeps the original as prof_sld-save.untilled.txt.
import sys, os, shutil
import numpy as np

MIX = {'g2', 'amfe3', 'gt', 'gr', 'sid', 'fes', 'cc', 'ka', 'gb', 'amsi', 'urea', 'fe2ad'}
d = sys.argv[1]; depth = float(sys.argv[2]) if len(sys.argv) > 2 else 0.16
f = os.path.join(d, 'prof', 'prof_sld-save.txt')
lines = open(f).read().splitlines()
head = lines[0]
names = ['z'] + [head[4+5*i:9+5*i].strip() for i in range((len(head.rstrip())-4)//5)]
x = np.array([[float(v) for v in l.split()] for l in lines[1:] if l.strip()])
if names[-1] != 'time' or len(names) != x.shape[1]:
    sys.exit(f'till.py: cannot parse header of {f}: {names}')
top = x[:, 0] < depth
for j, n in enumerate(names):
    if n in MIX: x[top, j] = x[top, j].mean()
shutil.copy(f, f.replace('.txt', '.untilled.txt'))
with open(f, 'w') as o:
    o.write(head+'\n')
    for r in x: o.write(''.join('%27.16E' % v for v in r)+'\n')
