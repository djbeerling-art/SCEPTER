# Collect global archetype runs (run_global.sh) into per-archetype, per-cell and per-nation season GHG tables.
# Usage: python3 paddy/global/collect.py <outdir> <cells.csv> [doses, default 10,20]   (from the repository root,
#   SCEPTER_OUT set as for the runs). Writes <outdir>/arch_ghg.csv, cell_ghg.csv, nation_ghg.csv.
# Per-ha values come from paddy/ghg_budget.py against each archetype's no-basalt run (one season, per ha);
# cell and nation totals multiply by harvested rice area (one season per year; multiple cropping not yet counted).
import sys, os, subprocess
import pandas as pd

D, cells = sys.argv[1], pd.read_csv(sys.argv[2])
doses = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else '10,20').split(',')]
O = os.path.join(os.environ['SCEPTER_OUT'], 'paddy', 'flood_')
COLS = ['ch4', 'n2o', 'co2deg', 'alk', 'calcite', 'diss', 'cations', 'd_ch4', 'd_n2o', 'cdr', 'net', 'net_pot']

rows = []
for a in pd.read_csv(os.path.join(D, 'archetypes.csv')).itertuples():
    runs = [f'{O}g_{a.arch}_{d}' for d in [0]+doses]
    if not all(os.path.isdir(r+'/flx') for r in runs): print('missing', a.arch); continue
    txt = subprocess.run(['python3', 'paddy/ghg_budget.py']+runs+['--dose='+','.join(map(str, doses))],
                         capture_output=True, text=True).stdout
    for l in [l.split() for l in txt.splitlines() if l.startswith('flood_')]:
        rows.append(dict(arch=a.arch, dose=float(l[1]), **dict(zip(COLS, map(float, l[2:])))))
ag = pd.DataFrame(rows); ag.to_csv(os.path.join(D, 'arch_ghg.csv'), index=False)

cm = pd.read_csv(os.path.join(D, 'cell_map.csv')).merge(ag, on='arch')
per_ha = COLS[:5]+COLS[6:]
for c in per_ha: cm[c+'_tot'] = cm[c]*cm.area          # per-ha x ha (kg or t, as in ghg_budget)
cell = cm.groupby(['cell', 'dose'])[[c+'_tot' for c in per_ha]+['area']].sum().reset_index()
cell = cell.merge(cells[['cell', 'lon', 'lat', 'iso3']], on='cell'); cell.to_csv(os.path.join(D, 'cell_ghg.csv'), index=False)
nat = cell.groupby(['iso3', 'dose'])[[c+'_tot' for c in per_ha]+['area']].sum().reset_index()
for c in ['d_ch4', 'd_n2o', 'cdr', 'net', 'net_pot']: nat[c+'_Mt'] = nat[c+'_tot']/1e6   # t CO2e -> Mt CO2e
nat.to_csv(os.path.join(D, 'nation_ghg.csv'), index=False)
print(nat[nat.dose > 0].groupby('dose')[['area']+[c+'_Mt' for c in ['d_ch4', 'd_n2o', 'cdr', 'net', 'net_pot']]].sum())
