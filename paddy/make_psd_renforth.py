# Writes a SCEPTER psdrain.in (sum of log-normal number PSDs) for rock dust with the particle size distribution of
# Renforth (2012, Int. J. Greenh. Gas Control 10, 229-243, Table 1; P80 ~10 um), with all sizes scaled by one factor
# so that the volume-weighted P80 (diameter, as SCEPTER's calc_p80 reports it) is the target.
# Usage: python3 paddy/make_psd_renforth.py [p80_um=100] [out=data/psdrain_renforth_p80-100um.in]
import sys, numpy as np
d_mm = np.array([0.0977,0.0691,0.0488,0.0345,0.0244,0.0173,0.0122,0.00863,0.004,0.002,0.001,0.0005,0.0001,0.00009])
ret = np.array([0.8,0.9,0.8,1.1,4.2,4.1,10.4,6.2,12.8,11.1,10.2,10.8,13.9,12.7])  # % retained (mass)
# each class holds material finer than its sieve size and coarser than the next: represent it at the geometric mean
d_lo = np.append(d_mm[1:], d_mm[-1]**2/d_mm[-2])
d_rep = np.sqrt(d_mm*d_lo)*1e-3  # m (diameter)
sig = 0.1  # log10 width of each component (SCEPTER bin width is ~0.12)
ps = np.linspace(np.log10(10e-9), np.log10(10e-3), 50)  # SCEPTER radius grid (log10 m)
def p80(scale):
    r = d_rep*scale/2; w = ret/100/r**3  # number weights giving the mass fractions
    psd = sum(wi*np.exp(-0.5*((ps-np.log10(ri))/sig)**2) for wi, ri in zip(w, r))
    cum = np.cumsum(psd*(10**ps)**3); cum /= cum[-1]
    dm = 2*10**ps*1e6
    i = np.where((cum[:-1] <= 0.8) & (cum[1:] >= 0.8))[0][-1]
    return dm[i] + (0.8-cum[i])*(dm[i+1]-dm[i])/(cum[i+1]-cum[i]), r, w
target = float(sys.argv[1]) if len(sys.argv) > 1 else 100.
out = sys.argv[2] if len(sys.argv) > 2 else f'data/psdrain_renforth_p80-{target:g}um.in'
lo, hi = 1., 100.
for _ in range(60):
    mid = np.sqrt(lo*hi)
    lo, hi = (mid, hi) if p80(mid)[0] < target else (lo, mid)
p, r, w = p80(mid)
with open(out, 'w') as f:
    f.write(f'** mean radius [m], standard deviation in log10 [-], weight [-]: Renforth 2012 Table 1 PSD x{mid:.3f}, P80 {p:.1f} um\n')
    for ri, wi in zip(r, w): f.write(f'{ri:.4e}\t{sig}\t{wi/w.max():.4e}\n')
print(f'scale x{mid:.3f}  P80 {p:.1f} um  ->  {out}')
