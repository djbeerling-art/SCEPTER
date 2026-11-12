"""Figures for the 27-soil arsenic response map (input: as_sweep_summary.csv from as_sweep_summary.py).
1. response_map_y<Y>.png: log10 competitive As(III) uptake ratio (Ki 300 uM), basalt / control, on a 3x3 grid of
   ferrihydrite x goethite, one panel per irrigation alkalinity; diverging scale centred on 0 (blue = uptake down,
   red = up), cells outlined where the class is 'uncertain (uptake parameters)'.
2. decomposition_y<Y>.png: log10 R_up = log10 R_As,eff - log10 T_eff per soil (As(III) mobilisation vs Si inhibition).
3. years.png: log10 R_up in the first vs last chain year per soil (1:1 line).
4. predictor.png: f/q capacity rule (uncorrected and pre-registered corrections A = 0.84x, B = 0.72x) vs the model's
   R_As, with log RMSE in the legend; R_inv vs f.
Usage: python3 paddy/arsenic/as_sweep_figures.py as_sweep_summary.csv [out dir, default .]"""
import sys, os, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
rows = list(csv.DictReader(open(sys.argv[1]))); out = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(out, exist_ok=True)
F = lambda r, k: float(r[k]) if r.get(k, '') not in ('', None) else np.nan
FE = sorted({float(r['fe3_wt']) for r in rows}); GT = sorted({float(r['gt_wt']) for r in rows})
AL = sorted({float(r['irralk']) for r in rows}); YR = sorted({int(r['year']) for r in rows})
INK, MUTED = '#1f2328', '#6e7781'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED})

def response_map(y):
    sel = [r for r in rows if int(r['year']) == y]
    v = np.array([np.log10(F(r, 'R_up_Ki300')) for r in sel]); lim = max(0.3, np.nanmax(np.abs(v)))
    fig, axs = plt.subplots(1, len(AL), figsize=(3.2*len(AL), 3.2), constrained_layout=True, squeeze=False)
    norm = TwoSlopeNorm(0, -lim, lim)
    for ax, a in zip(axs[0], AL):
        Z = np.full((len(GT), len(FE)), np.nan); U = np.zeros_like(Z, bool)
        for r in sel:
            if float(r['irralk']) != a: continue
            i, j = GT.index(float(r['gt_wt'])), FE.index(float(r['fe3_wt']))
            Z[i, j] = np.log10(F(r, 'R_up_Ki300')); U[i, j] = r['class'].startswith('uncertain')
        cm = plt.get_cmap('RdBu_r').copy(); cm.set_bad('#d0d7de')
        im = ax.imshow(np.ma.masked_invalid(Z), origin='lower', cmap=cm, norm=norm)
        for i in range(len(GT)):
            for j in range(len(FE)):
                if np.isfinite(Z[i, j]):
                    ax.text(j, i, f'{10**Z[i, j]:.2f}', ha='center', va='center', color=INK if abs(Z[i, j]) < 0.6*lim else 'white', fontsize=8)
                else: ax.text(j, i, 'pending', ha='center', va='center', color=MUTED, fontsize=7)
                if U[i, j]: ax.add_patch(plt.Rectangle((j-.5, i-.5), 1, 1, fill=False, ec=INK, lw=1.5, ls='--'))
        ax.set_xticks(range(len(FE)), [f'{100*x:g}' for x in FE]); ax.set_yticks(range(len(GT)), [f'{100*x:g}' for x in GT])
        ax.set_xlabel('ferrihydrite (wt%)'); ax.set_ylabel('goethite (wt%)'); ax.set_title(f'irrigation alkalinity x{a:g}', color=INK)
    cb = fig.colorbar(im, ax=axs[0].tolist(), shrink=0.8); cb.set_label('log10 As(III) uptake ratio (basalt/control)')
    fig.suptitle(f'Year {y}: blue = uptake suppressed, red = raised; R = 1 at white; grey = run pending; dashed = uptake-parameter uncertainty', fontsize=9, color=INK)
    fig.savefig(os.path.join(out, f'response_map_y{y}.png'), dpi=200); plt.close(fig)

def decomposition(y):
    sel = sorted([r for r in rows if int(r['year']) == y], key=lambda r: F(r, 'log10R_up'))
    x = np.arange(len(sel)); fig, ax = plt.subplots(figsize=(9, 3.4), constrained_layout=True)
    ax.bar(x-0.2, [F(r, 'log10R_As_eff') for r in sel], 0.4, color='#c0504d', label='log10 R_As(III) (mobilisation)')
    ax.bar(x+0.2, [-F(r, 'log10T_eff') for r in sel], 0.4, color='#4f81bd', label='- log10 T (Si inhibition)')
    ax.plot(x, [F(r, 'log10R_up') for r in sel], 'o', color=INK, ms=4, label='log10 R_uptake (sum)')
    ax.axhline(0, color=MUTED, lw=0.8); ax.set_xticks(x, [r['soil'].replace('as_', '') for r in sel], rotation=90, fontsize=6)
    ax.set_ylabel('log10 ratio'); ax.legend(frameon=False, fontsize=7); ax.set_title(f'Year {y}: log R_uptake = log R_As(III) - log T', fontsize=9)
    fig.savefig(os.path.join(out, f'decomposition_y{y}.png'), dpi=200); plt.close(fig)

def years():
    if len(YR) < 2: return
    a = {r['soil']: F(r, 'log10R_up') for r in rows if int(r['year']) == YR[0]}
    b = {r['soil']: F(r, 'log10R_up') for r in rows if int(r['year']) == YR[-1]}
    k = [s for s in a if s in b]; fig, ax = plt.subplots(figsize=(3.6, 3.6), constrained_layout=True)
    ax.plot([a[s] for s in k], [b[s] for s in k], 'o', color='#4f81bd', ms=5)
    lo, hi = min(min(a.values()), min(b.values()), -0.1), max(max(a.values()), max(b.values()), 0.1)
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8); ax.axhline(0, color=MUTED, lw=0.5, ls=':'); ax.axvline(0, color=MUTED, lw=0.5, ls=':')
    ax.set_xlabel(f'log10 R_uptake, year {YR[0]}'); ax.set_ylabel(f'log10 R_uptake, year {YR[-1]}')
    fig.savefig(os.path.join(out, 'years.png'), dpi=200); plt.close(fig)

def predictor():
    m = np.array([F(r, 'R_As') for r in rows]); p = np.array([F(r, 'R_As_fq') for r in rows]); f = np.array([F(r, 'f_fh') for r in rows])
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.5), constrained_layout=True)
    ax = axs[0]
    for fac, lab, col in ((1, 'f/q rule', '#9e9e9e'), (0.84, 'A: 0.84 x rule', '#4f81bd'), (0.72, 'B: 0.72 x rule', '#c0504d')):
        e = np.sqrt(np.nanmean((np.log10(fac*p) - np.log10(m))**2))
        ax.plot(fac*p, m, 'o', ms=4, color=col, label=f'{lab} (log RMSE {e:.3f})')
    lo, hi = np.nanmin([m, 0.72*p]), np.nanmax([m, p]); ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('predicted R_As(III)'); ax.set_ylabel('model R_As(III)'); ax.legend(frameon=False, fontsize=7)
    ax = axs[1]; ax.plot(f, [F(r, 'R_inv') for r in rows], 'o', ms=4, color='#4f81bd')
    ax.axhline(0.84, color=MUTED, lw=0.8, ls='--'); ax.set_xlabel('ferrihydrite share of As(III) capacity, f'); ax.set_ylabel('As(III) inventory ratio R_inv')
    fig.savefig(os.path.join(out, 'predictor.png'), dpi=200); plt.close(fig)

for y in YR: response_map(y); decomposition(y)
years(); predictor()
print('figures in', out)
