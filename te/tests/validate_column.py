"""Benchmark the one-way coupled TE column (chem/column.py) against PHREEQC
TRANSPORT.

An acid soil column (20 cells; exchanger with Al + Hfo + gibbsite; Ni, Cu, Pb
and Cd loaded) is flushed
with EW-like Ca-Mg-HCO3 water carrying no metals. PHREEQC solves the full
problem. The Python column takes only the host state (pH, water, major-element
system totals, Al included) from PHREEQC and moves and partitions the
metals itself. Pass if
total, dissolved and exchangeable metal agree with PHREEQC everywhere to
1% of the column's initial per-cell inventory (dissolved: of the initial
dissolved), and the effluent to 1% of its peak.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from chem import thermo, system, solver, column, phreeqc_ref as pr  # noqa: E402

OUT = os.path.join(HERE, 'out')
MAJ = ['Ca', 'Mg', 'K', 'Na', 'Cl', 'C', 'Al']
METALS = ['Ni', 'Cu', 'Pb', 'Cd']
NCELL, SHIFTS = 20, 60
HFO = dict(weak=0.025, strong=6e-4, mass=11.0, area=600.0)
SPEC = dict(
    ncell=NCELL, shifts=SHIFTS,
    inflow={'pH': 7.5, 'Ca': 1.5e-3, 'Mg': 5e-4, 'C(4)': 4e-3, 'O(0)': 5e-4,
            'charge_on': 'C(4)'},
    soil={'pH': 5.0, 'Ca': 3e-4, 'Mg': 1e-4, 'K': 5e-5, 'Na': 2e-4, 'Cl': 1e-3,
          'C(4)': 1e-4, 'O(0)': 5e-4, 'Al': 1e-4, 'Ni': 2e-6, 'Cu': 1e-6,
          'Pb': 2e-7, 'Cd': 5e-8, 'charge_on': 'Cl'},
    phases=[('Gibbsite', 0, 0.01)],
    cec=0.1, hfo=HFO, elements=MAJ + METALS, metals=METALS,
    ex_species={m: [m + 'X2'] for m in METALS})


def host_state(r):
    return {'pH': r['pH'], 'W': r['W'],
            'totals': {e: r['SYS_' + e] for e in MAJ}}


def main():
    t0 = time.time()
    init, rows = pr.run_transport(SPEC)
    t_ph = time.time() - t0
    db = thermo.read([pr.DB, pr.EXTRA])
    sy = system.build(db, MAJ + METALS, valences={'C': ['C(4)'], 'Cu': ['Cu(2)']},
                      exchangers=['X'], surfaces=['Hfo_w', 'Hfo_s'], phases=['Gibbsite'])
    sv = solver.Speciator(sy)
    host = [[host_state(init)] * NCELL]
    host += [[host_state(rows[(k, i + 1)]) for i in range(NCELL)]
             for k in range(1, SHIFTS + 1)]
    init_tot = [{m: init['SYS_' + m] for m in METALS} for _ in range(NCELL)]
    surf = {'Hfo': solver.Surface({'Hfo_w': HFO['weak'], 'Hfo_s': HFO['strong']},
                                  HFO['area'], HFO['mass'])}
    t0 = time.time()
    res = column.run(sv, host, METALS, init_tot, {}, cec={'X': SPEC['cec']},
                     surfaces=surf, minerals=['Gibbsite'])
    t_py = time.time() - t0

    lines = ['# TE column vs PHREEQC TRANSPORT', '',
             f'{NCELL} cells, {SHIFTS} shifts, pure advection; exchanger '
             f"{SPEC['cec']} mol + Hfo ({HFO['mass']} g) + gibbsite per cell; acid soil "
             'solution (pH 5) flushed with Ca-Mg-HCO3 water (pH 7.5), no metal inflow.',
             f'PHREEQC {t_ph:.1f} s, Python column {t_py:.1f} s.', '',
             '| metal | pool | max abs err / scale | scale (mol) |', '|---|---|---|---|']
    worst = 0.0
    for j, m in enumerate(METALS):
        ph = {p: np.zeros((SHIFTS, NCELL)) for p in ('total', 'dissolved', 'exchange')}
        for k in range(1, SHIFTS + 1):
            for i in range(NCELL):
                r = rows[(k, i + 1)]
                ph['total'][k - 1, i] = r['SYS_' + m]
                ph['dissolved'][k - 1, i] = r['TOT_' + m] * r['W']
                ph['exchange'][k - 1, i] = r['EX_' + m] * r['W']
        for pool, scale in (('total', init['SYS_' + m]),
                            ('exchange', init['SYS_' + m]),
                            ('dissolved', init['TOT_' + m] * init['W'])):
            py = res.pools[pool][1:, :, j]
            e = np.max(np.abs(py - ph[pool])) / scale
            worst = max(worst, e)
            lines.append(f'| {m} | {pool} | {e:.2e} | {scale:.3g} |')
        eff_ph = ph['dissolved'][:, -1]
        eff_py = res.pools['dissolved'][1:, -1, j]
        e = np.max(np.abs(eff_py - eff_ph)) / max(eff_ph.max(), 1e-30)
        worst = max(worst, e)
        lines.append(f'| {m} | effluent (last cell) | {e:.2e} | {eff_ph.max():.3g} |')
        lines.append(f'| {m} | mass left after {SHIFTS} shifts | '
                     f'PHREEQC {ph["total"][-1].sum()/init["SYS_"+m]/NCELL:.3f}, '
                     f'Python {res.pools["total"][-1, :, j].sum()/init["SYS_"+m]/NCELL:.3f} '
                     f'of initial | |')
    ok = worst < 0.01
    lines += ['', f'Worst error {worst:.2e}: {"PASS" if ok else "FAIL"}']
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, 'column_validation.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    _plot(rows, res, init)
    return 0 if ok else 1


def _plot(rows, res, init):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(1, 1 + len(METALS), figsize=(4 * (1 + len(METALS)), 3.6))
    x = np.arange(1, NCELL + 1)
    for k, ls in ((10, '-'), (30, '--'), (60, ':')):
        ax[0].plot(x, [rows[(k, i)]['pH'] for i in x], 'k' + ls, label=f'shift {k}')
        for j, (m, c) in enumerate(zip(METALS, ('C0', 'C3', 'C2', 'C4'))):
            ax[1 + j].plot(x, [rows[(k, i)]['SYS_' + m] / init['SYS_' + m] for i in x],
                           color=c, ls=ls, label=f'PHREEQC, shift {k}')
            ax[1 + j].plot(x, res.pools['total'][k, :, j] / init['SYS_' + m], 'o',
                           color=c, ms=3, mfc='none')
    ax[0].set(xlabel='cell', ylabel='pH (host)')
    ax[0].legend(frameon=False, fontsize=8)
    for j, m in enumerate(METALS):
        ax[1 + j].set(xlabel='cell', ylabel=f'total {m} / initial', title=f'{m} (o = Python)')
    ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'column_profiles.png'), dpi=130)


if __name__ == '__main__':
    sys.exit(main())
