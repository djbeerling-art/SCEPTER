# Season greenhouse-gas budget of basalt-amended paddy runs against their no-basalt control.
# Usage: python3 paddy/ghg_budget.py <control run dir> <run dir> [<run dir> ...] [--days=113] [--dose=<t/ha>,...]
# Per run (season totals, per ha):
#   CH4 emission  = diffusion + ebullition + unoxidised half of plant transport (as in calib_metrics.py)
#   N2O emission  = diffusion + plant transport
#   CO2 degassing = net pco2 exchange across the top (positive = to the atmosphere)
#   ALK export    = alkalinity leached through the base: the EW CO2 removal term, counted 1 mol CO2 per eq
#                   (assumes the alkalinity is kept downstream)
#   calcite       = net calcite precipitation in the column (1 mol CO2 stored per mol)
# Net GHG change vs control (t CO2-eq/ha) = dCH4*27 + dN2O*273 - dALK export - dcalcite (AR6 GWP100, non-fossil CH4).
import sys,numpy as np

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt  = dict(a[2:].split('=') for a in sys.argv[1:] if a.startswith('--') and '=' in a)
days = float(opt.get('days',113.))
tend = days/365.
BASALT = ['fo','ab','an','fa','dp','kfs']

def integ(run,name):
    f = f'{run}/flx/flx_{name}.txt'
    h = open(f).readline().split()
    d = np.loadtxt(f,skiprows=1)
    t = np.concatenate([[0],d[:,0]])
    dt = np.diff(np.minimum(t,tend))
    return {h[i]:np.sum(d[:,i]*dt) for i in range(1,len(h))}

def season(run):
    c = integ(run,'gas-pch4')
    ch4 = c['dif'] + c['ch4eb'] + 0.5*c['ch4pl']
    n = integ(run,'gas-pn2o')
    o = {'CH4 (kg/ha)':ch4*16.04*10,'N2O (kg N/ha)':(n['dif']+n.get('n2opl',0))*28.01*10}
    o['CO2 degassing (t/ha)'] = integ(run,'gas-pco2')['dif']*0.4401
    o['ALK export (t CO2/ha)'] = integ(run,'co2sp-ALK')['adv']*0.4401
    try: o['calcite (t CO2/ha)'] = -integ(run,'sld-cc')['cc']*0.4401   # cc column: dissolution positive
    except OSError: o['calcite (t CO2/ha)'] = 0.
    diss = 0.
    for m in BASALT:
        try: diss += integ(run,f'sld-{m}')[m]
        except (OSError,KeyError): pass
    o['basalt dissolved (mol/m2)'] = diss
    return o

ctrl = season(args[0])
doses = [float(x) for x in opt['dose'].split(',')] if 'dose' in opt else [np.nan]*len(args[1:])
keys = list(ctrl)
print('run'.ljust(14)+'dose'.rjust(6)+''.join(k.split(' (')[0][:10].rjust(12) for k in keys)+'  dCO2e_CH4  dCO2e_N2O  CDR  net')
for run,dose in zip([args[0]]+args[1:],[0.]+doses):
    o = season(run)
    dch4 = (o['CH4 (kg/ha)']-ctrl['CH4 (kg/ha)'])/1e3*27
    dn2o = (o['N2O (kg N/ha)']-ctrl['N2O (kg N/ha)'])*44.01/28.01/1e3*273
    cdr  = (o['ALK export (t CO2/ha)']-ctrl['ALK export (t CO2/ha)']) + (o['calcite (t CO2/ha)']-ctrl['calcite (t CO2/ha)'])
    print(run.split('/')[-1][:14].ljust(14)+f'{dose:6.0f}'+''.join(f'{o[k]:12.4g}' for k in keys)
          +f'{dch4:11.3f}{dn2o:11.3f}{cdr:6.3f}{dch4+dn2o-cdr:7.3f}')
