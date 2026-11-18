"""Summarise the arsenic sweep (paddy/arsenic/run_as_sweep.sh): one row per soil and treatment year.

Diagnostics (basalt B vs control C, root zone 0-ZROOT m, root-density weighted exp(-z/0.1) as ch4pl/siup):
  R_As      season-mean porewater As(III) B/C                         (geochemical mobilisation)
  R_Si, dpH season-mean Si B/C and pH B-C                             (dissolved Si response, alkalinisation)
  T         (1 + Si_B/Ki)/(1 + Si_C/Ki) from season means, Ki = 300 uM (maximum mobilisation compatible with suppression)
  R_up_Ki*  competitive uptake ratio from TIME-RESOLVED values: integral over time and roots of As3/(1 + Si/Ki), B/C,
            for Ki = 150, 300, 600 uM (linear in As(III): Km_As >> porewater As(III))
  log10 decomposition: log R_up(Ki 300) = log R_As,eff - log T_eff, where R_As,eff = time-integrated As(III) B/C and
            T_eff = R_As,eff / R_up (the Si-inhibition term from the time-resolved series)
  R_up_sim  As(III) uptake from the model's as3up flux, B/C (int_flx_aq-as3: time-averaged rate x time)
  budget    per run (mol As m-2, As(III)+As(V)): percolating irrigation input (rain.in concs x q x t), net transport
            gain (-(adv+dif) x t: percolation input + floodwater exchange (irrigation-renewed floodwater, lateral
            discharge) - drainage at the base), plant uptake, change in sorbed and dissolved inventory, closure residual
  f_fh, q   ferrihydrite-capacity predictor (equilibrium, linear sorption, conserved exchangeable As(III)):
            Kd_fh, Kd_gt = sorbed/dissolved As(III) on ferrihydrite / goethite sites (as3sp equation, root-weighted
            time means); f_fh = Kd_fh/(Kd_fh+Kd_gt) in the control; q = Kd_fh,C/Kd_fh,B (Si and pH), q_Si = Si part only;
            R_As_fq = 1/((1-f_fh) + f_fh/q); R_As_cap = (1+Kd_C)/(1+Kd_B) (also carries the goethite pH term)
  R_inv     exchangeable As(III) inventory (dissolved + sorbed) B/C; at local equilibrium R_As = R_inv x R_cap, so
            R_inv carries what the predictor leaves out (As(V) reduction, As(III) oxidation, transport, uptake)
  class     'suppress' / 'mobilise' if every uptake estimate (Ki 150/300/600 and simulated) is below / above 1;
            'uncertain (uptake parameters)' if they disagree: uptake-parameter uncertainty only, not Fe sorption,
            As(V) reduction or grain partitioning
  grain iAs not modelled yet (translocation, Lsi1/Lsi2 regulation, node sequestration) -> empty
Usage: python3 paddy/arsenic/as_sweep_summary.py <SCEPTER_OUT>/paddy [dose | si] > as_sweep_summary.csv
       ('si': silica-only test of run_si_test.sh, control _sc vs silica _ss)
"""
import sys, os, glob, re, csv
import numpy as np
_trap = getattr(np, "trapezoid", None) or getattr(np, "trapz")

root = sys.argv[1]; dose = sys.argv[2] if len(sys.argv) > 2 else '20'
ARMC, ARMB = ('sc', 'ss') if dose == 'si' else ('d0', f'd{dose}')
ZROOT, ZR = 0.2, 0.1
KIS = (150.0, 300.0, 600.0)

A3 = (5.1882, -0.0224, -0.0171, 3.8028, 0.0104, 0.8685)   # as3sp constants (kd_constants.json)

def profs(run, kind):
    """(header, array) per output time. Some files have fused header names (fixed-width fields run together): such a
    file takes the header of another file of the same kind whose names match its columns, else its own line split at
    the field width that gives the right number of names."""
    out = []
    for f in sorted(glob.glob(os.path.join(run, 'prof', f'{kind}-[0-9]*.txt'))):
        raw = open(f).readline().rstrip('\n')
        try: a = np.atleast_2d(np.loadtxt(f, skiprows=1))
        except ValueError: continue
        if a.size == 0: continue
        out.append([raw.split(), a, raw])
    good = [o[0] for o in out if len(o[0]) == o[1].shape[1]]
    for o in out:
        if len(o[0]) == o[1].shape[1]: continue
        nc = o[1].shape[1]; cand = [h for h in good if len(h) == nc]
        if cand: o[0] = cand[0]; continue
        st = o[2].index('z') + 2                         # fields start after 'z ' (till.py: 4, width 5)
        for wd in (5, 6, 7, 8):
            h = ['z'] + [o[2][st+wd*i:st+wd*(i+1)].strip() for i in range((len(o[2].rstrip()) - st + wd - 1)//wd)]
            if len(h) == nc and h[-1] == 'time': o[0] = h; break
    return [(h, a) for h, a, _ in out if len(h) == a.shape[1]]

def sites(run):
    try: lines = open(os.path.join(run, 'as_sites.in')).read().splitlines()
    except OSError: lines = []
    v = []
    for l in lines:   # header line (text) then 'f [fgt]'; keep the first all-numeric line
        try: v = [float(x) for x in l.split()]; break
        except ValueError: continue
    return (v[0] if v else 1.0), (v[1] if len(v) > 1 else 0.02)

def kd3(run):
    """As(III) sorption constants of the run: ./as_kd.in (written by --askd, e.g. a high-phosphate refit) or built-in"""
    try: return tuple(float(x) for x in open(os.path.join(run, 'as_kd.in')).read().splitlines()[1].split()[:6])
    except (OSError, IndexError, ValueError): return A3

def series(run):
    """time-resolved As(III), Si (M), pH and As(III) Kd on ferrihydrite and goethite sites over the root zone; Kd and
    inventory are NaN where the solid profiles were pruned (earlier chain years keep porewater only, KEEPPH)"""
    aq, sl, bs = profs(run, 'prof_aq'), profs(run, 'prof_sld'), profs(run, 'bsd')
    if len(aq) < 3: return None
    full = len(sl) == len(aq) == len(bs)
    fs, fg = sites(run); a, b, c, aS, bS, n = kd3(run)
    z = aq[0][1][:, 0]; w = np.exp(-z/ZR)*(z <= ZROOT); w = w/w.sum()
    o = {k: [] for k in ('t', 'as3', 'si', 'ph', 'kfh', 'kgt', 'dsi', 'inv3')}
    for i, (h, x) in enumerate(aq):
        g = lambda hh, xx, k: xx[:, hh.index(k)]
        ph, si = g(h, x, 'ph'), g(h, x, 'si'); pc, sc = np.clip(ph, 5, 8) - 7, np.clip(si, 2e-5, 2e-3)
        K = 10**(a + b*pc + c*pc*pc); D = 1 + 10**(aS + bS*pc)*sc**n
        o['t'].append(x[0, h.index('time')]); o['as3'].append(g(h, x, 'as3')); o['si'].append(si); o['ph'].append(ph); o['dsi'].append(D)
        if full:
            (hs, s), (hb, bb) = sl[i], bs[i]
            th = 1e3*g(hb, bb, 'poro')*g(hb, bb, 'sat')
            o['kfh'].append(K*fs*0.2*g(hs, s, 'amfe3')/D/th); o['kgt'].append(K*fg*g(hs, s, 'gt')/th if 'gt' in hs else 0*K)
            o['inv3'].append(g(h, x, 'as3')*th + g(hs, s, 'as3ad'))   # exchangeable As(III), mol m-3
        else:
            for k in ('kfh', 'kgt', 'inv3'): o[k].append(np.full_like(ph, np.nan))
    o = {k: np.array(v) for k, v in o.items()}; o['w'] = w
    return o

def tint(s, x):
    """time integral (trapezoid) of a root-weighted profile quantity x[t, z]"""
    v = (x*s['w']).sum(axis=1)
    return _trap(v, s["t"]) if len(v) > 1 else v[0]

def mean(s, x): return tint(s, x)/(s['t'][-1] - s['t'][0]) if len(s['t']) > 1 else (x*s['w']).sum()

def intflx(run, name):
    f = os.path.join(run, 'flx', f'int_flx_{name}.txt')
    if not os.path.exists(f): return None
    h = open(f).readline().split()
    try: a = np.atleast_2d(np.loadtxt(f, skiprows=1))
    except ValueError: return None
    if a.size == 0 or a.shape[1] != len(h): return None
    return {k: a[-1, i]*a[-1, 0] for i, k in enumerate(h) if k != 'time'}   # time-averaged rate x time = total

def budget(run):
    a3, a5, s3, s5 = (intflx(run, n) for n in ('aq-as3', 'aq-as5', 'sld-as3ad', 'sld-as5ad'))
    if a3 is None or a5 is None: return {}
    conc = {l.split()[0]: float(l.split()[1]) for l in open(os.path.join(run, 'rain.in')).readlines()[1:] if len(l.split()) > 1}
    fr = open(os.path.join(run, 'frame.in')).readlines()
    ttot, q = float(fr[3].split()[0]), float(fr[15].split()[0])
    inp = (conc.get('as3', 0) + conc.get('as5', 0))*1e3*q*ttot
    gain = -(a3.get('adv', 0) + a3.get('dif', 0) + a5.get('adv', 0) + a5.get('dif', 0))
    upt = a3.get('as3up', 0)
    dsorb = (s3['tflx'] + s5['tflx']) if s3 and s5 else np.nan
    ddis = a3['tflx'] + a5['tflx']
    return dict(input_perc=inp, net_transport_gain=gain, uptake=upt, d_sorbed=dsorb, d_dissolved=ddis,
                residual=gain - upt - dsorb - ddis)

w = csv.writer(sys.stdout)
cols = ['soil', 'fe3_wt', 'gt_wt', 'irralk', 'year', 'pH_C', 'dpH', 'Si_C_uM', 'R_Si', 'As3_C_uM', 'R_As', 'T',
        'R_up_Ki150', 'R_up_Ki300', 'R_up_Ki600', 'R_up_sim', 'f_fh', 'q', 'q_Si', 'R_As_fq', 'R_As_cap', 'R_inv', 'log10R_As_eff', 'log10T_eff', 'log10R_up', 'class',
        'grain_iAs_ratio']
for tr in ('C', 'B'):
    cols += [f'{tr}_{k}_molm2' for k in ('input_perc', 'net_transport_gain', 'uptake', 'd_sorbed', 'd_dissolved', 'residual')]
w.writerow(cols)
for d0 in sorted(glob.glob(os.path.join(root, f'flood_as_f*_{ARMC}_y*s'))):
    m = re.search(rf'flood_(as_f([\d.]+)_g([\d.]+)_a([\d.]+))_{ARMC}_y(\d+)s$', d0)
    if not m: continue
    tag, f, g, a, y = m.groups()
    db = os.path.join(root, f'flood_{tag}_{ARMB}_y{y}s')
    if not os.path.isdir(db): continue
    c, b = series(d0), series(db)
    if c is None or b is None:
        print(f'# skipped {os.path.basename(d0)} / {os.path.basename(db)}: incomplete porewater profiles', file=sys.stderr); continue
    sc, sb = mean(c, c['si']), mean(b, b['si'])
    ra = mean(b, b['as3'])/mean(c, c['as3']); T = (1 + sb*1e6/300)/(1 + sc*1e6/300)
    rup = {ki: tint(b, b['as3']/(1 + b['si']*1e6/ki))/tint(c, c['as3']/(1 + c['si']*1e6/ki)) for ki in KIS}
    ra_eff = tint(b, b['as3'])/tint(c, c['as3']); t_eff = ra_eff/rup[300.0]
    kfc, kgc, kfb, kgb = mean(c, c['kfh']), mean(c, c['kgt']), mean(b, b['kfh']), mean(b, b['kgt'])
    f_fh = kfc/(kfc + kgc); q = kfc/kfb; q_si = mean(b, b['dsi'])/mean(c, c['dsi'])
    r_fq = 1/((1 - f_fh) + f_fh/q); r_cap = (1 + kfc + kgc)/(1 + kfb + kgb)
    bc, bb = budget(d0), budget(db)
    nr = lambda v, k=3: '' if not np.isfinite(v) else round(v, k)
    us = bb.get('uptake', np.nan)/bc['uptake'] if bc.get('uptake') else np.nan
    est = list(rup.values()) + ([us] if np.isfinite(us) else [])
    cls = 'suppress' if max(est) < 1 else ('mobilise' if min(est) > 1 else 'uncertain (uptake parameters)')
    row = [tag, f, g, a, y, round(mean(c, c['ph']), 2), round(mean(b, b['ph']) - mean(c, c['ph']), 2), round(sc*1e6),
           round(sb/sc, 2), '%.4g' % (mean(c, c['as3'])*1e6), round(ra, 3), round(T, 3),
           *[round(rup[k], 3) for k in KIS], round(us, 3) if np.isfinite(us) else '',
           nr(f_fh), nr(q), round(q_si, 3), nr(r_fq), nr(r_cap), nr(mean(b, b['inv3'])/mean(c, c['inv3'])),
           round(np.log10(ra_eff), 3), round(np.log10(t_eff), 3), round(np.log10(rup[300.0]), 3), cls, '']
    for bd in (bc, bb):
        row += ['%.4g' % bd.get(k, np.nan) for k in ('input_perc', 'net_transport_gain', 'uptake', 'd_sorbed', 'd_dissolved', 'residual')]
    w.writerow(row)
