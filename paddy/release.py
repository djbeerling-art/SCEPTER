# Cation release by primary-mineral dissolution over a season (mmol m-2 d-1), Uchibayashi-style.
# Usage: python3 release.py <run dir> [species ...]; primary minerals = rock minerals present in the flx header
import sys, numpy as np
d=sys.argv[1]; sp=sys.argv[2:] or ['mg','na','ca','k','si']
PRIM=['by','la','an','ab','agt','aug','dps','en','fo','fa','mgn','ilm','olg','and','kfs','sm','cpx','opx','hb','qtz']
for s in sp:
    try: f=open(f'{d}/flx/flx_aq-{s}.txt')
    except OSError: continue
    h=f.readline().split(); a=np.loadtxt(f)
    t=a[:,0]; cols=[i for i,n in enumerate(h) if n in PRIM]
    r=a[:,cols].sum(axis=1)                    # mol m-2 yr-1 (SCEPTER sign: dissolution source)
    tot=np.trapezoid(r,t)                          # mol m-2
    print(f'{s:3s} {",".join(h[i] for i in cols):25s} {tot*1e3/(t[-1]*365.25):9.2f} mmol/m2/d over {t[-1]*365.25:.0f} d')
