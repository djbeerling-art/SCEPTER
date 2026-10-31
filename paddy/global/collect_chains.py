# Collect global multi-year archetype chains (stanage_handoff/run_arch_v22.sh, tags g22_<arch>_<dose>_y<N>{s,f})
# into per-archetype, per-cell and per-nation year-round GHG / CDR tables, using the same budget as the GMD v22
# archetype figure (v22/figcdr22.py): season + fallow each year, treatment minus the year-matched dose-0 chain,
# calcite counted as CDR only up to the rock cations not exported as alkalinity (cumulative cap), GWP100 (CH4 27, N2O 273).
# Usage (repository root):
#   python3 paddy/global/collect_chains.py <run dir = $SCEPTER_OUT/paddy> <archetypes.csv> <cell_map.csv> <cells.csv> <outdir>
#       [--doses=10,20] [--years=10] [--tag=g22] [--alk=gross|pda (default pda)]
# Writes <outdir>/arch_years.csv (per archetype, dose, year), arch_cum.csv (cumulative), cell_cum.csv, nation_cum.csv.
# Per-ha values are t CO2e/ha (cumulative over the chain); totals multiply by harvested rice area (one season per year).
import sys, os
import numpy as np, pandas as pd

pos = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
R, ARCH, CMAP, CELLS, OUT = pos
doses = [int(x) for x in opt.get('doses', '10,20').split(',')]
NY, TAG = int(opt.get('years', 10)), opt.get('tag', 'g22')
sys.argv = ['x', '--days=500']                       # whole fallow window for the fallow half of each year
exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0])   # defines season()
AK = 'ALK export (t CO2/ha)' if opt.get('alk') == 'gross' else 'potentially durable ALK (t CO2/ha)'
K = lambda n: n*44/28*0.273                          # kg N2O-N/ha -> t CO2e/ha
os.makedirs(OUT, exist_ok=True)

yrs, cum = [], []
for a in pd.read_csv(ARCH).itertuples():
    run = lambda d, y, x: f'{R}/flood_{TAG}_{a.arch}_{d}_y{y}{x}'
    if not os.path.isdir(run(0, NY, 'f')): print('missing control', a.arch); continue
    ctrl = {y: {x: season(run(0, y, x)) for x in 'sf'} for y in range(1, NY+1)}
    for d in doses:
        if not os.path.isdir(run(d, NY, 'f')): print('missing', a.arch, d); continue
        t = dict(ch4s=0., ch4f=0., n2o=0., alk=0., cc=0., cat=0.)
        for y in range(1, NY+1):
            b, c = {x: season(run(d, y, x)) for x in 'sf'}, ctrl[y]
            r = dict(arch=a.arch, dose=d, year=y,
                     ch4_ctrl=sum(c[x]['CH4 (kg/ha)'] for x in 'sf'), ch4_trt=sum(b[x]['CH4 (kg/ha)'] for x in 'sf'),
                     n2o_ctrl=sum(c[x]['N2O (kg N/ha)'] for x in 'sf'), n2o_trt=sum(b[x]['N2O (kg N/ha)'] for x in 'sf'),
                     ch4s=(b['s']['CH4 (kg/ha)']-c['s']['CH4 (kg/ha)'])*27e-3,
                     ch4f=(b['f']['CH4 (kg/ha)']-c['f']['CH4 (kg/ha)'])*27e-3,
                     n2o=K(sum(b[x]['N2O (kg N/ha)']-c[x]['N2O (kg N/ha)'] for x in 'sf')),
                     alk=sum(b[x][AK]-c[x][AK] for x in 'sf'),
                     cc=sum(b[x]['calcite (t CO2/ha)']-c[x]['calcite (t CO2/ha)'] for x in 'sf'),
                     cat=sum(b[x]['cation release (t CO2/ha)'] for x in 'sf'))
            yrs.append(r)
            for k in t: t[k] += r[k]
        t['cc'] = min(t['cc'], max(t['cat']-t['alk'], 0.)); t['ret'] = t['cat']-t['alk']-t['cc']
        t['net'] = t['ch4s']+t['ch4f']+t['n2o']-t['alk']-t['cc']
        t['net_pot'] = t['ch4s']+t['ch4f']+t['n2o']-t['cat']
        cum.append(dict(arch=a.arch, dose=d, years=NY, **t))
pd.DataFrame(yrs).to_csv(os.path.join(OUT, 'arch_years.csv'), index=False)
ag = pd.DataFrame(cum); ag.to_csv(os.path.join(OUT, 'arch_cum.csv'), index=False)
if ag.empty: sys.exit('no complete archetypes')

V = ['ch4s', 'ch4f', 'n2o', 'alk', 'cc', 'cat', 'ret', 'net', 'net_pot']
cm = pd.read_csv(CMAP).merge(ag, on='arch')
for v in V: cm[v+'_t'] = cm[v]*cm.area                # t CO2e/ha x ha = t CO2e over the chain
cell = cm.groupby(['cell', 'dose'])[[v+'_t' for v in V]+['area']].sum().reset_index()
cell = cell.merge(pd.read_csv(CELLS, usecols=['cell', 'lon', 'lat', 'iso3']), on='cell')
cell.to_csv(os.path.join(OUT, 'cell_cum.csv'), index=False)
nat = cell.groupby(['iso3', 'dose'])[[v+'_t' for v in V]+['area']].sum().reset_index()
for v in V: nat[v+'_Mt'] = nat[v+'_t']/1e6
nat.to_csv(os.path.join(OUT, 'nation_cum.csv'), index=False)
tot = nat.groupby('dose')[['area']+[v+'_Mt' for v in V]].sum()
tot['area'] /= 1e6
print(f'Cumulative {NY}-yr totals (area Mha; Mt CO2e, negative = removal):'); print(tot.round(1).to_string())
