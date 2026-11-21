# Season metrics of a paddy column run, for comparison with field targets.
# Usage: python3 paddy/calib_metrics.py <run directory> [season length, days; default 113 (IPCC 2019 global default)]
#        [fertiliser N, kg N/ha]
# CH4 emission = diffusion out of the top + the half of plant-transported CH4 that is not oxidised in the
# rhizosphere (ch4pl stoichiometry) + ebullition.
import sys,glob,numpy as np

run  = sys.argv[1]
days = float(sys.argv[2]) if len(sys.argv) > 2 else 113.
tend = days/365.

def flx(name):
    f = f'{run}/flx/flx_{name}.txt'
    h = open(f).readline().split()
    d = np.loadtxt(f,skiprows=1)
    t = np.concatenate([[0],d[:,0]])
    dt = np.diff(np.minimum(t,tend))
    return {h[i]:np.sum(d[:,i]*dt) for i in range(1,len(h))}, d, h

c,d,h = flx('gas-pch4')
plant = c.get('ch4pl',0)*0.5
ebul  = c.get('ch4eb',0)
dif   = c.get('dif',0)
emis  = dif + plant + ebul
out = {
    'CH4 production (mol/m2)': -c.get('g2ch4',0),
    'CH4 emission (g/m2/season)': emis*16.04,
    'CH4 emission (kg/ha/d)': emis*16.04*10/days,
    'plant share of emission': plant/emis if emis > 0 else np.nan,
    'ebullition share': ebul/emis if emis > 0 else np.nan,
}
# emission rate time series: onset = first time the rate reaches 10% of its in-season maximum
rate = 0.5*d[:,h.index('ch4pl')] + d[:,h.index('ch4eb')] + d[:,h.index('dif')]
ins = d[:,0] <= tend
if rate[ins].max() > 0:
    out['CH4 onset (d, 10% of max rate)'] = d[ins,0][np.argmax(rate[ins] >= 0.1*rate[ins].max())]*365

# porewater maxima over the season (profiles)
def prof_max(kind,col):
    best,tbest = 0.,np.nan
    for f in sorted(glob.glob(f'{run}/prof/prof_{kind}-*.txt')):
        hh = open(f).readline().split()
        if col not in hh: return np.nan,np.nan
        dd = np.loadtxt(f,skiprows=1)
        if dd[0,-1] > tend: continue
        v = dd[:,hh.index(col)].max()
        if v > best: best,tbest = v,dd[0,-1]
    return best,tbest*365
fe2,tfe2 = prof_max('aq','fe2')
ch4,tch4 = prof_max('gas','pch4')
out['max porewater Fe2+ (mM)'] = fe2*1e3
out['  day of Fe2+ max'] = tfe2
out['max porewater CH4 (mM, kH 1.4e-3 M/atm)'] = ch4*1.4e-3*1e3
try:
    cn,_,_ = flx('gas-pn2o')
    # losses are positive in the flux files: dif = exchange across the top, n2opl = plant transport
    n2o = cn.get('dif',0) + cn.get('n2opl',0)
    out['N2O emission (mg N/m2/season)'] = n2o*28e3
    out['N2O plant share'] = cn.get('n2opl',0)/n2o if n2o > 0 else np.nan
    if len(sys.argv) > 3:  # fertiliser N [kg N/ha] for the emission factor (no background subtraction)
        out['N2O-N / fertiliser N (%)'] = n2o*28e3/(float(sys.argv[3])*100)*100
except OSError:
    pass
for k,v in out.items(): print(f'{k:42s} {v:.4g}')
