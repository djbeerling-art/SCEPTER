# Morris elementary-effects screen of the controls on the season CH4 (and net GHG) response to basalt.
# design:  python3 paddy/screen/morris.py design <dir> [--r=8] -> <dir>/design.csv, jobs_spin.txt, jobs_ctrl.txt, jobs_dose.txt
#          (one job line = tag + driver options; identical spin-ups and controls are shared between design points)
# analyse: python3 paddy/screen/morris.py analyse <dir>  -> <dir>/morris_<metric>.csv (mu*, sigma per factor)
# Factors and ranges (Morris p=4 levels): soil reducible Fe(III), OM input, soil pH (through irrigation alkalinity,
# 0 -> pH ~5.0, 2 -> ~6.7 after spin-up), soil CEC, temperature, rock reducible Fe(III), basalt (5 SCEPTER rocks),
# dose, water management (continuous flooding with mid-season drain / AWD), green-rust bioavailability (g2gr krxn2,
# the least constrained calibration constant).
import sys, os, hashlib
import numpy as np, pandas as pd

FACT = [('soilfe3', 0.004, 0.03), ('om', 250, 800), ('irralk', 0, 2), ('cec', 8, 30), ('temp', 20, 30),
        ('basfe3', 0, 0.05), ('rock', 0, 4), ('dose', 5, 30), ('awd', 0, 1), ('grbio', 0.02, 0.3)]
ROCKS = ['Hillhouse', 'Oregon', 'BlueRidge', 'Tichum', 'basalt']
SPIN = ['soilfe3', 'om', 'irralk', 'cec', 'temp']
CTRL = SPIN+['awd', 'grbio']
PROB = dict(num_vars=len(FACT), names=[f[0] for f in FACT], bounds=[[f[1], f[2]] for f in FACT])
KRXN = 'paddy/krxn_calibrated_v4.in'

def key(row, cols): return hashlib.md5(','.join(f'{row[c]:.6g}' for c in cols).encode()).hexdigest()[:8]

def opts(r):
    o = f'--fe3={r.soilfe3:.4g} --om={r.om:.0f} --irralk={r.irralk:.3g} --cec={r.cec:.3g} --temp={r.temp:.3g}'
    return o

def design(d, nr):
    from SALib.sample import morris
    X = morris.sample(PROB, N=nr, num_levels=4, seed=1)
    df = pd.DataFrame(X, columns=PROB['names'])
    df['rock'] = df.rock.round().astype(int); df['awd'] = df.awd.round().astype(int)
    df['spin'] = [key(r, SPIN) for _, r in df.iterrows()]
    df['ctrl'] = [key(r, CTRL) for _, r in df.iterrows()]
    df['tag'] = [f'm_{i:03d}' for i in range(len(df))]
    df.to_csv(os.path.join(d, 'design.csv'), index=False)
    os.makedirs(os.path.join(d, 'krxn'), exist_ok=True)
    base = open(KRXN).read().splitlines()
    for g in sorted(set(df.grbio)):
        with open(os.path.join(d, 'krxn', f'gr{g:.4f}.in'), 'w') as f:
            f.write('\n'.join(f'g2gr -999 {g:.4f}' if l.startswith('g2gr') else l for l in base)+'\n')
    sp, ct, ds = {}, {}, []
    for r in df.itertuples():
        kf = os.path.join(d, 'krxn', f'gr{r.grbio:.4f}.in'); w = '--awd' if r.awd else '--phd=14'
        sp[f'ms_{r.spin}'] = opts(r)
        ct[f'mc_{r.ctrl}'] = f'--rstrt=ms_{r.spin} --k={kf} {w} --basalt=0 {opts(r)}'
        ds.append(f'{r.tag} --rstrt=ms_{r.spin} --k={kf} {w} --basalt={r.dose:.3g} --rock={ROCKS[r.rock]} '
                  f'--basfe3={r.basfe3:.4g} {opts(r)}')
    for n, lines in (('spin', [f'{k} {v}' for k, v in sp.items()]), ('ctrl', [f'{k} {v}' for k, v in ct.items()]),
                     ('dose', ds)):
        open(os.path.join(d, f'jobs_{n}.txt'), 'w').write('\n'.join(lines)+'\n')
    print(f'{len(df)} design points: {len(sp)} spin-ups, {len(ct)} controls, {len(ds)} basalt runs')

def analyse(d):
    import subprocess
    from SALib.analyze import morris
    df = pd.read_csv(os.path.join(d, 'design.csv'))
    O = os.path.join(os.environ['SCEPTER_OUT'], 'paddy', 'flood_')
    res = []
    for r in df.itertuples():
        txt = subprocess.run(['python3', 'paddy/ghg_budget.py', O+'mc_'+r.ctrl, O+r.tag, f'--dose={r.dose:.3g}'],
                             capture_output=True, text=True).stdout
        l = [l.split() for l in txt.splitlines() if l.startswith('flood_')]
        if len(l) < 2: res.append([np.nan]*5); continue
        c, b = [float(x) for x in l[0][2:]], [float(x) for x in l[1][2:]]
        res.append([100*(b[0]/c[0]-1), b[7], b[8], b[9], b[10]])   # % CH4, dCH4 CO2e, dN2O CO2e, CDR, net
    out = pd.DataFrame(res, columns=['ch4_pct', 'd_ch4', 'd_n2o', 'cdr', 'net'])
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
