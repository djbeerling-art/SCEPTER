# Season greenhouse-gas budget of basalt-amended paddy runs against their no-basalt control.
# Usage: python3 paddy/ghg_budget.py <control run dir> <run dir> [<run dir> ...] [--days=113] [--dose=<t/ha>,...]
# Per run (season totals, per ha):
#   CH4 emission  = diffusion + ebullition + unoxidised half of plant transport (as in calib_metrics.py)
#   N2O emission  = diffusion + plant transport
#   CO2 degassing = net pco2 exchange across the top (positive = to the atmosphere)
#   ALK export    = alkalinity leached through the base: the EW CO2 removal term, counted 1 mol CO2 per eq
#                   (assumes the alkalinity is kept downstream)
#   calcite       = net calcite precipitation in the column (1 mol CO2 stored per mol)
#   cation release = base-cation charge released by basalt dissolution: the CDR potential, realised as the cations
#                   leave the soil (in one season most stay on the exchange complex)
# Net GHG change vs control (t CO2-eq/ha, AR6 GWP100, non-fossil CH4):
#   realised  = dCH4*27 + dN2O*273 - dALK export - dcalcite, with dcalcite capped at the basalt's cation release
#               less dALK export (calcite beyond that is not rock-derived CDR)
#   potential = dCH4*27 + dN2O*273 - cation release
import sys,numpy as np

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt  = dict(a[2:].split('=') for a in sys.argv[1:] if a.startswith('--') and '=' in a)
days = float(opt.get('days',113.))
tend = days/365.
warned = set()

def integ(run,name):
    f = f'{run}/flx/flx_{name}.txt'
    h = open(f).readline().split()
    d = np.loadtxt(f,skiprows=1)
    t = np.concatenate([[0],d[:,0]])
    if d[-1,0] < tend*0.999 and run not in warned:   # the model can stop early with a normal exit status
        warned.add(run); print(f'WARNING: {run} ends at {d[-1,0]*365:.1f} d, before the {days:.0f}-d season',file=sys.stderr)
    dt = np.diff(np.minimum(t,tend))
    return {h[i]:np.sum(d[:,i]*dt) for i in range(1,len(h))}

def season(run):
    c = integ(run,'gas-pch4')
    ch4 = c['dif'] + c.get('ch4eb',0) + 0.5*c.get('ch4pl',0)
    n = integ(run,'gas-pn2o')
    o = {'CH4 (kg/ha)':ch4*16.04*10,'N2O (kg N/ha)':(n['dif']+n.get('n2opl',0))*28.01*10}
    o['CO2 degassing (t/ha)'] = integ(run,'gas-pco2')['dif']*0.4401
    o['ALK export (t CO2/ha)'] = integ(run,'co2sp-ALK')['adv']*0.4401
    try: o['calcite (t CO2/ha)'] = -integ(run,'sld-cc')['cc']*0.4401   # cc column: dissolution positive
    except OSError: o['calcite (t CO2/ha)'] = 0.
    # rock minerals = those in the run's dust.in (calcite, Fe(III) oxide and inert excluded); their base-cation charge
    # release is read from the aqueous flux files (2 Ca + 2 Mg + Na + K), so any rock recipe works
    rock = [l.split()[0] for l in open(f'{run}/dust.in').readlines()[1:] if l.split() and l.split()[0] not in ('cc','amfe3','inrt')]
    diss = eq = 0.
    for m in rock:
        try: diss += integ(run,f'sld-{m}')[m]
        except (OSError,KeyError): pass
    for sp,z in [('ca',2),('mg',2),('na',1),('k',1)]:
        try: c = integ(run,f'aq-{sp}')
        except OSError: continue
        eq -= z*sum(c.get(m,0.) for m in rock)   # mineral columns: dissolution supplies the solute (negative)
    o['basalt dissolved (mol/m2)'] = diss
    o['cation release (t CO2/ha)'] = eq*0.4401   # CDR potential once the cations leave the soil
    return o

ctrl = season(args[0])
doses = [float(x) for x in opt['dose'].split(',')] if 'dose' in opt else [np.nan]*len(args[1:])
keys = list(ctrl)
print('run'.ljust(14)+'dose'.rjust(6)+''.join(k.split(' (')[0][:10].rjust(12) for k in keys)+'  dCO2e_CH4  dCO2e_N2O  CDR  net  net_pot')
for run,dose in zip([args[0]]+args[1:],[0.]+doses):
    o = season(run)
    dch4 = (o['CH4 (kg/ha)']-ctrl['CH4 (kg/ha)'])/1e3*27
    dn2o = (o['N2O (kg N/ha)']-ctrl['N2O (kg N/ha)'])*44.01/28.01/1e3*273
    dalk = o['ALK export (t CO2/ha)']-ctrl['ALK export (t CO2/ha)']
    dcc  = o['calcite (t CO2/ha)']-ctrl['calcite (t CO2/ha)']
    # extra calcite counts as CDR only up to the base-cation charge the basalt released (less what already left as
    # alkalinity): beyond that, its Ca and alkalinity come from the exchange complex and Fe/SO4 reduction, which is
    # not rock-derived and reverses on drainage (found in warm, Fe-rich archetypes, where a small pH lift tips the
    # soil into calcite precipitation of several t CO2/ha)
    cdr  = dalk + min(dcc, max(o['cation release (t CO2/ha)']-dalk, 0.)) if run != args[0] else 0.
    print(run.split('/')[-1][:14].ljust(14)+f'{dose:6.0f}'+''.join(f'{o[k]:12.4g}' for k in keys)
          +f' {dch4:10.3f} {dn2o:10.3f} {cdr:6.3f} {dch4+dn2o-cdr:7.3f} {dch4+dn2o-o["cation release (t CO2/ha)"]:8.3f}')
