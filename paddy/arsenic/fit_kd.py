"""SCEPTER Step 1 sorption: trace As on ferrihydrite weak sites as a site-normalised distribution coefficient
with Si competition, fitted separately for As(III) and As(V) to the PHREEQC DLM grid (grid_phreeqc.csv).

  sorbed_i (mol/L pw) = S * c_i * K_i / (1 + KSi_i * cSi^n_i),   S = 0.2 * Fe(III) as ferrihydrite (mol/L pw)
  log K_i   = a_i + b_i (pH-7) + c_i (pH-7)^2
  log KSi_i = aS_i + bS_i (pH-7)
As is trace relative to sites (<5% loading in the grid), so As(III) and As(V) do not compete with each other.
Mass balance for one species: c_i = T_i / (1 + S K_i / (1 + KSi_i cSi^n_i))."""
import csv, json, sys, numpy as np
from scipy.optimize import least_squares
GRID = sys.argv[1] if len(sys.argv) > 1 else 'grid_phreeqc.csv'
OUTJ = sys.argv[2] if len(sys.argv) > 2 else 'kd_constants.json'
rows = list(csv.DictReader(open(GRID)))
D = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
d = D['pH'] - 7; si = D['si_mM']*1e-3; S = 0.2*D['fe_mol']
T = {'as3': D['as_mM']*1e-3*(1-D['as5_frac']), 'as5': D['as_mM']*1e-3*D['as5_frac']}
obs = {'as3': D['as3_uM']*1e-6, 'as5': D['as5_uM']*1e-6}

def pred(p, sp):
    a, b, c, aS, bS, n = p
    K = 10**(a + b*d + c*d*d); KS = 10**(aS + bS*d)
    return T[sp]/(1 + S*K/(1 + KS*si**n))

out = {}
for sp, p0 in (('as3', [3, 0, 0, 4, 0, 0.8]), ('as5', [6, -0.5, 0, 4, 0, 0.8])):
    f = least_squares(lambda p: np.log10(pred(p, sp)/obs[sp]), p0)
    r = np.log10(pred(f.x, sp)/obs[sp])
    out[sp] = dict(zip(['a','b','c','aS','bS','n'], f.x.round(4).tolist()))
    out[sp]['rmse_log10'] = round(float(np.sqrt(np.mean(r**2))), 3)
    out[sp]['max_abs_log10'] = round(float(abs(r).max()), 3)
    print(sp, out[sp])
json.dump(out, open(OUTJ, 'w'), indent=1)

# behaviour checks against the prototype scans (base case: Fe 0.233, As 0.208 mM, 5% As(V))
def c_of(sp, pH, si_mM, fe=0.233):
    p = [out[sp][k] for k in ['a','b','c','aS','bS','n']]
    a, b, c, aS, bS, n = p; dd = pH - 7
    K = 10**(a + b*dd + c*dd*dd); KS = 10**(aS + bS*dd)
    return 1/(1 + 0.2*fe*K/(1 + KS*(si_mM*1e-3)**n))  # dissolved fraction
print('As(III) dissolved ratio, Si 0.3 -> 1.5 mM, pH 6.8: %.2f' % (c_of('as3', 6.8, 1.5)/c_of('as3', 6.8, 0.3)))
print('As(III) dissolved ratio, 80%% ferrihydrite loss: %.2f' % (c_of('as3', 6.8, 0.3, 0.233*0.2)/c_of('as3', 6.8, 0.3)))
print('As(III) dissolved ratio, 90%% ferrihydrite loss: %.2f' % (c_of('as3', 6.8, 0.3, 0.233*0.1)/c_of('as3', 6.8, 0.3)))
print('As(V) dissolved ratio, pH 6 -> 7.5 (Si 0.3): %.2f' % (c_of('as5', 7.5, 0.3)/c_of('as5', 6.0, 0.3)))
