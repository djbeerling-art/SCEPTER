#!/usr/bin/env python3
# Alkalinity and reduced-species fluxes across depth planes, to see how much of the alkalinity leaving the puddled
# layer survives transport through the plough pan and an oxidising subsoil.
# Usage: python3 paddy/alk_flux.py <run dir> [<run dir> ...] [--q=0.5] [--z=0.2,0.24,bottom] [--ts=<csv>] [--rhob=1350]
#        [--fe3in=<mol m-2 yr-1, as given to the driver>] [--zsupp=0.3]
#   (several run dirs, e.g. the seasons and fallows of a chain, are concatenated in the order given)
#
# Planes are cell interfaces (0.2 = base of the puddled layer, 0.24 = base of the pan); the advective flux through a
# plane is q times the concentration of the cell above it (upwind). At each plane the advective flux q*c (mol m-2) is integrated over the profile records; diffusion is neglected (it
# is small next to percolation and is cut further in the pan). Concentrations are dissolved totals.
#   A_carb  carbonate alkalinity from charge balance: 2Ca + 2Mg + Na + K + 2Fe2+ + NH4+ - 2SO4 - Cl - NO3 - HS-
#           (Al and silicate anions neglected; at the column base this can be checked against flx_co2sp-ALK)
#   Fe2+    total dissolved Fe(II);  NH4+  total ammonia N;  S_T  total dissolved sulfide;  HS-  free bisulfide
#   PDA     potentially durable alkalinity = A_carb - 2 Fe2+ - 2 NH4+ - (2 S_T - HS-)
#           i.e. the carbonate alkalinity left once the exported reduced species are oxidised: Fe2+ to Fe(OH)3
#           (2 H+ per Fe, assuming oxidation, hydrolysis and precipitation), NH4+ by complete nitrification (2 H+ per
#           N; denitrification or uptake downstream would return some), sulfide to sulfate (2 eq per S, less the HS-
#           that is already alkalinity). A diagnostic, not verified permanent CO2 removal.
# Fe check (planes 2 and 3): Fe2+ flux retained between them vs the change in Fe held in that layer: solids (amfe3,
# gt, siderite, FeS, green rust), exchangeable Fe(II) (cmol/kg x bulk density --rhob, kg/m3) and porewater Fe(II).
import sys, os, glob
import numpy as np, pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--'))
q = float(opt.get('q', 0.5))
planes = opt.get('z', '0.2,0.24,bottom').split(',')

def records(run, kind):
    fs = sorted(glob.glob(os.path.join(run, 'prof', f'{glob.escape(kind)}-[0-9][0-9][0-9].txt')))
    fs = [f for f in fs if not f.endswith('-000.txt')]   # record 000 is the initial state, with a different header
    return [pd.read_csv(f, sep=r'\s+') for f in fs]

T, DT, TOT, FREE, GAS, SLD, ADS = [], [], [], [], [], [], []
for run in args:
    tot, free, gas = records(run, 'prof_aq(tot)'), records(run, 'prof_aq'), records(run, 'prof_gas(aq)')
    if not tot:
        sys.exit(f'no prof_aq(tot) records in {run} (rerun with KEEPALK=1)')
    n = min(len(tot), len(free), len(gas))
    t = np.array([r['time'].iloc[0] for r in tot[:n]])   # each run's clock starts at 0: chain the runs end to end
    DT += list(np.diff(np.concatenate([[0.], t]))); T += list(t + (T[-1] if T else 0.)); TOT += tot[:n]; FREE += free[:n]; GAS += gas[:n]
    SLD.append(records(run, 'prof_sld')); ADS.append(records(run, 'prof_aq(ads)'))
T = np.array(T); z = TOT[0]['z'].values
dt = np.array(DT)   # rectangle rule on the records

def g(r, k): return r[k].values*1e3 if k in r else np.zeros(len(r))   # mol/L -> mol per m3 of water
fe = np.array([g(r, 'fe2') for r in TOT])
nh4 = np.array([g(r, 'pnh3') for r in GAS])
st = np.array([g(r, 'hs') for r in TOT])
hs = np.array([g(r, 'hs') for r in FREE])
acarb = np.array([2*g(r, 'ca') + 2*g(r, 'mg') + g(r, 'na') + g(r, 'k') - 2*g(r, 'so4') - g(r, 'cl') - g(r, 'no3')
                  for r in TOT]) + 2*fe + nh4 - hs
pda = acarb - 2*fe - 2*nh4 - (2*st - hs)
series = {'A_carb': acarb, 'Fe2': fe, 'NH4': nh4, 'S_T': st, 'HS': hs, 'PDA': pda}

print(f'{len(args)} run(s), {len(T)} records, t {T[0]:.3f}-{T[-1]:.3f} yr, q {q} m/yr; cumulative fluxes, mol m-2')
print('plane_m ' + ''.join(f'{k:>10}' for k in series) + '  PDA_tCO2/ha')
ts = {'time': T}; cum = {}
for p in planes:
    zp = z[-1] + 1e-6 if p == 'bottom' else float(p)
    iu = np.where(z < zp)[0][-1]   # upwind cell: percolation is downward, so the flux through a plane carries the cell above
    row = {}
    for k, c in series.items():
        f = q*c[:, iu]   # mol m-2 yr-1 through the plane
        ts[f'{k}@{zp:g}'] = f; row[k] = np.sum(f*dt)
    cum[zp] = row
    print(f'{min(zp, z[-1]):7.3f} ' + ''.join(f'{row[k]:10.4f}' for k in series) + f'  {row["PDA"]*0.4401:10.4f}')

zs = sorted(cum)
if len(zs) >= 3 and SLD[0] and SLD[-1]:
    a, b = zs[1], zs[2]
    def inv(r, cols):
        m = (r['z'] > a) & (r['z'] < b); dz = np.gradient(r['z'].values)[m]
        return sum(w*float(np.sum(r[c].values[m]*dz)) for c, w in cols if c in r)
    fe3c = [('amfe3', 1), ('gt', 1)]
    fetc = fe3c + [('sid', 1), ('fes', 1), ('gr', 3)]   # green rust here holds 1 Fe(II) + 2 Fe(III)
    rhob = float(opt.get('rhob', 1350.))
    def exch(r):   # exchangeable Fe(II), cmol/kg -> mol per m3 bulk
        m = (r['z'] > a) & (r['z'] < b); dz = np.gradient(r['z'].values)[m]
        return float(np.sum(r['fe2'].values[m]*1e-2*rhob*dz))
    def held(i):   # i = 0 first record, -1 last
        s_, a_ = SLD[i][i], (ADS[i][i] if ADS[i] else None)
        return inv(s_, fetc), (exch(a_) if a_ is not None else np.nan)
    (s0, e0), (s1, e1) = held(0), held(-1)
    lost = cum[a]['Fe2'] - cum[b]['Fe2']
    zs_ = float(opt.get('zsupp', 0.3))   # Fe(III) supplied with the OM input (--fe3in) is spread as exp(-z/zsupp)/zsupp
    fin = float(opt.get('fe3in', 0)) * (T[-1] - T[0]) * (np.exp(-a/zs_) - np.exp(-b/zs_))
    print(f'Fe check {a:g}-{b:g} m: Fe2+ flux retained {lost:.4f} + Fe(III) input with OM {fin:.4f} = {lost + fin:.4f} mol m-2; change in solid Fe {s1 - s0:.4f} '
          f'(Fe(III) oxides {inv(SLD[-1][-1], fe3c) - inv(SLD[0][0], fe3c):.4f}), exchangeable Fe(II) {e1 - e0:.4f}, '
          f'sum {s1 - s0 + e1 - e0:.4f} (porewater storage not counted)')
if 'ts' in opt:
    pd.DataFrame(ts).to_csv(opt['ts'], index=False); print('time series ->', opt['ts'])
