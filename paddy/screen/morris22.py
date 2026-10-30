# Morris elementary-effects screen of the controls on the season CH4 (and net GHG) response to basalt.
# design:  python3 paddy/screen/morris22.py design <dir> [--r=8] -> <dir>/design.csv, jobs_spin.txt, jobs_ctrl.txt, jobs_dose.txt
#          (one job line = tag + driver options; identical spin-ups and controls are shared between design points)
# analyse: python3 paddy/screen/morris22.py analyse <dir>  -> <dir>/morris_<metric>.csv (mu*, sigma per factor)
# Factors and ranges (Morris p=4 levels): soil reducible Fe(III), OM input, soil pH (through irrigation alkalinity,
# 0 -> pH ~5.0, 2 -> ~6.7 after spin-up), soil CEC, temperature, rock reducible Fe(III), basalt (5 SCEPTER rocks),
# dose, water management (continuous flooding with pre-harvest drain / AWD), strength w of the Si feedback on root O2
# loss (frozen at 0.2; least constrained constant, --rtsi) and percolation q (0.25-1.1 m/yr = 0.7-3 mm/d).
# v22 port (morris22.py): frozen v22 physics in the plough-pan column; each point = 3-yr drained spin-up, 7-yr settle
# (control years, own water management), then one year (season + fallow) of control and of basalt applied once.
import sys, os, hashlib
import numpy as np, pandas as pd

FACT = [('soilfe3', 0.004, 0.03), ('om', 250, 900), ('irralk', 0, 2), ('cec', 8, 30), ('temp', 20, 32),
        ('basfe3', 0, 0.05), ('rock', 0, 4), ('dose', 5, 30), ('awd', 0, 1), ('siw', 0.1, 0.4), ('q', 0.25, 1.1)]
ROCKS = ['Hillhouse', 'Oregon', 'BlueRidge', 'Tichum', 'basalt']
SPIN = ['soilfe3', 'om', 'irralk', 'cec', 'temp', 'q']
CTRL = SPIN+['awd', 'siw']
PROB = dict(num_vars=len(FACT), names=[f[0] for f in FACT], bounds=[[f[1], f[2]] for f in FACT])
SETL = SPIN+['awd']

def key(row, cols): return hashlib.md5(','.join(f'{row[c]:.6g}' for c in cols).encode()).hexdigest()[:8]

def opts(r):
    o = f'--fe3={r.soilfe3:.4g} --om={r.om:.0f} --irralk={r.irralk:.3g} --cec={r.cec:.3g} --temp={r.temp:.3g} --q={r.q:.3g}'
    return o

def design(d, nr):
    from SALib.sample import morris
    X = morris.sample(PROB, N=nr, num_levels=4, seed=1)
    df = pd.DataFrame(X, columns=PROB['names'])
    df['rock'] = df.rock.round().astype(int); df['awd'] = df.awd.round().astype(int)
    df['spin'] = [key(r, SPIN) for _, r in df.iterrows()]
    df['setl'] = [key(r, SETL) for _, r in df.iterrows()]
    df['ctrl'] = [key(r, CTRL) for _, r in df.iterrows()]
    df['tag'] = [f'm22_{i:03d}' for i in range(len(df))]
    df.to_csv(os.path.join(d, 'design.csv'), index=False)
    sp, st, ct, ds = {}, {}, {}, []
    for r in df.itertuples():
        w = '--awd ' if r.awd else ''; si = f'--rtsi=2e-3:2.7e-4:{r.siw:.3g}'
        sp[f'ms22_{r.spin}'] = opts(r)
        st[f'mv22_{r.setl}'] = f'ms22_{r.spin} {w}{opts(r)}'
        ct[f'mc22_{r.ctrl}'] = f'mv22_{r.setl} 0 {w}{si} {opts(r)}'
        ds.append(f'{r.tag} mv22_{r.setl} {r.dose:.3g} {w}{si} --rock={ROCKS[r.rock]} --basfe3={r.basfe3:.4g} {opts(r)}')
    for n, lines in (('spin', [f'{k} {v}' for k, v in sp.items()]), ('setl', [f'{k} {v}' for k, v in st.items()]),
                     ('ctrl', [f'{k} {v}' for k, v in ct.items()]), ('dose', ds)):
        open(os.path.join(d, f'jobs_{n}.txt'), 'w').write('\n'.join(lines)+'\n')
    print(f'{len(df)} design points: {len(sp)} spin-ups, {len(st)} settles, {len(ct)} controls, {len(ds)} basalt runs')

def analyse(d):
    # metrics per design point, one year (season + fallow), basalt run minus its year-matched control:
    #   base_ch4 = control season CH4 (kg/ha): baseline emissions (depends only on soil, climate, water and w)
    #   ch4_pct  = season CH4 change (%); d_ch4 = year CH4 change (t CO2e/ha, GWP100 27)
    #   net      = dCH4 + dN2O - dPDA - dcalcite (calcite capped at cations not exported), t CO2e/ha
    from SALib.analyze import morris
    sys.argv = ['x', '--days=500', '--alk=pda']
    g = {}; exec(open('paddy/ghg_budget.py').read().split('ctrl = season')[0], g); season = g['season']
    AK = 'potentially durable ALK (t CO2/ha)'
    df = pd.read_csv(os.path.join(d, 'design.csv'))
    O = os.path.join(os.environ['SCEPTER_OUT'], 'paddy', 'flood_')
    res = []
    for r in df.itertuples():
        try:
            a = {x: season(f'{O}mc22_{r.ctrl}_y1{x}') for x in 'sf'}; b = {x: season(f'{O}{r.tag}_y1{x}') for x in 'sf'}
        except (OSError, ValueError, IndexError):
            res.append([np.nan]*6); continue
        dd = lambda k: sum(b[x][k]-a[x][k] for x in 'sf')
        cat = sum(b[x]['cation release (t CO2/ha)'] for x in 'sf'); pda = dd(AK); cc = min(dd('calcite (t CO2/ha)'), max(cat-pda, 0))
        dch4 = dd('CH4 (kg/ha)')*27e-3; dn2o = dd('N2O (kg N/ha)')*44/28*0.273
        res.append([a['s']['CH4 (kg/ha)'], 100*(b['s']['CH4 (kg/ha)']/a['s']['CH4 (kg/ha)']-1), dch4, dn2o, pda+cc, dch4+dn2o-pda-cc])
    out = pd.DataFrame(res, columns=['base_ch4', 'ch4_pct', 'd_ch4', 'd_n2o', 'inorg_cdr', 'net'])
    pd.concat([df, out], axis=1).to_csv(os.path.join(d, 'results.csv'), index=False)
    X = df[PROB['names']].values.astype(float)
    for m in out:
        y = out[m].values
        if np.isnan(y).any(): print(m, 'has', np.isnan(y).sum(), 'failed runs; filled with the mean'); y = np.where(np.isnan(y), np.nanmean(y), y)
        s = morris.analyze(PROB, X, y, num_levels=4)
        t = pd.DataFrame({k: s[k] for k in ('names', 'mu', 'mu_star', 'sigma', 'mu_star_conf')}).sort_values('mu_star', ascending=False)
        t.to_csv(os.path.join(d, f'morris_{m}.csv'), index=False); print(m); print(t.round(3).to_string(index=False))

if __name__ == '__main__':
    d = sys.argv[2]; os.makedirs(d, exist_ok=True)
    nr = int(next((x.split('=')[1] for x in sys.argv if x.startswith('--r=')), 8))
    design(d, nr) if sys.argv[1] == 'design' else analyse(d)
