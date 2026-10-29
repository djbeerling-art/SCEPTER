# Annual diagnostics for a v22 10-yr chain vs its year-matched control: season CH4 production / oxidation / emission
# (kg CH4/ha), Fe(III) reduction (mol/m2), plough-layer (0-0.2 m) reactive Fe(III) inventory at the end of each year
# (start-state profile + cumulative layer budgets), and Ca/Mg/ALK fluxes through 0.2 m and the column base (mol/m2/yr).
# Usage: python3 annual22.py <outdir> <ctrl chain> <trt chain> <start state run> [<start state ctrl>]
import sys, glob, numpy as np, pandas as pd
O, C, T, S = sys.argv[1:5]; SC = sys.argv[5] if len(sys.argv) > 5 else S
def integ(d, f, cols):
    x = pd.read_csv(f'{d}/flx/{f}', sep=r'\s+'); t = np.concatenate([[0], x.time.values]); dt = np.diff(t)
    return {c: float((x[c].values*dt).sum()) if c in x else 0. for c in cols}
def zint(d, f, cols, zmax):
    try: x = pd.read_csv(f'{d}/flx/zint_{f}.txt', sep=r'\s+')
    except FileNotFoundError: return 0.
    return float(x[cols].sum(axis=1)[x.z < zmax].sum())
def fe3_0(run, zmax=0.2):
    f = sorted(glob.glob(f'{O}/flood_{run}/prof/prof_sld-[0-9]*.txt'))[-1]; p = pd.read_csv(f, sep=r'\s+')
    dz = np.diff(np.concatenate([[0], (p.z.values[:-1]+p.z.values[1:])/2, [p.z.values[-1]+(p.z.values[-1]-p.z.values[-2])/2]]))
    return float((p.amfe3*dz)[p.z < zmax].sum())
def ch4(d):
    g = integ(d, 'flx_gas-pch4.txt', ['g2ch4', 'ch4o2', 'ch4rx', 'dif', 'ch4pl', 'ch4eb']); k = 16e-3*1e4
    fe = integ(d, 'flx_aq-fe2.txt', ['g2fe3', 'g2gr', 'g2gt'])
    return -g['g2ch4']*k, (g['ch4o2']+g['ch4rx'])*k, (g['dif']+g['ch4pl']+g['ch4eb'])*k, -sum(fe.values())
def flux(d, sp, zmax): return zint(d, sp, ['adv', 'dif'], zmax)    # net transport loss (export) from the layer above zmax
inv = {C: fe3_0(SC), T: fe3_0(S)}
print(f'{"yr":>2s} | {"CH4 emis c/t (kg)":>18s} {"d%":>5s} | {"prod c/t":>13s} {"ox c/t":>13s} | {"FeRed c/t":>11s} | {"Fe(III) 0-20 c/t":>16s} | '
      f'{"dCa out 20cm/base":>14s} {"dMg 20cm/base":>14s} {"dALK 20cm/base":>15s}')
for y in range(1, 11):
    r = {}
    for n in (C, T):
        s, f = f'{O}/flood_{n}_y{y}s', f'{O}/flood_{n}_y{y}f'
        p, o, e, fr = ch4(s)
        inv[n] += sum(zint(d, 'sld-amfe3', ['tflx'], 0.2) for d in (s, f))
        fl = {sp: [sum(flux(d, sp, z) for d in (s, f)) for z in (0.2, 9)] for sp in ('aq-ca', 'aq-mg', 'co2sp-ALK')}
        r[n] = (p, o, e, fr, inv[n], fl)
    c, t = r[C], r[T]; d = {sp: [t[5][sp][i]-c[5][sp][i] for i in (0, 1)] for sp in t[5]}
    print(f'{y:2d} | {c[2]:8.1f} {t[2]:8.1f} {100*(t[2]/c[2]-1):+5.0f} | {c[0]:6.0f} {t[0]:6.0f} {c[1]:6.0f} {t[1]:6.0f} | {c[3]:5.2f} {t[3]:5.2f} | '
          f'{c[4]:7.1f} {t[4]:7.1f} | {d["aq-ca"][0]:+6.2f} {d["aq-ca"][1]:+6.2f} {d["aq-mg"][0]:+6.2f} {d["aq-mg"][1]:+6.2f} '
          f'{d["co2sp-ALK"][0]:+6.2f} {d["co2sp-ALK"][1]:+6.2f}')
