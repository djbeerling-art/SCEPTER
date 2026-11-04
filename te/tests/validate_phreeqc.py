"""Benchmark te/chem against PHREEQC (minteq.v4 + te_extra.dat).

Batch pH edges (pH 4 to 8) at fixed pe, each solved by PHREEQC and by the
Python speciator from the same system totals. Writes
te/tests/out/validation.csv, validation.md and pH_edges.png.

    python3 te/tests/validate_phreeqc.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from chem import thermo, system, solver, humic, phreeqc_ref as pr  # noqa: E402

OUT = os.path.join(HERE, 'out')
PH = [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
METALS = ['Ni', 'Cu', 'Zn', 'Cd', 'Pb']
HFO = dict(weak=2e-4, strong=5e-6, mass=0.089, area=600)   # 1 mmol Fe, D&M 1990

# HMO: Tonkin et al. (2004); 0.1 g (2.1 mmol sites/g, 64% X sites)
HMO = dict(x=0.1 * 2.1e-3 * 0.64, y=0.1 * 2.1e-3 * 0.36, mass=0.1, area=746)

db = thermo.read([pr.DB, pr.EXTRA])


def base(**kw):
    c = dict(T=25, W=1.0, pe=12, pco2_si=None, cec=None, hfo=HFO, minerals=[],
             solution={'Na': 0.01, 'Cl': 0.011}, pH=PH)
    c.update(kw)
    return c


CASES = []
for m in METALS:
    CASES.append(('edge_' + m, base(solution={'Na': 0.01, 'Cl': 0.011, m: 1e-6}),
                  {'metals': [m]}))
# Cr: PHREEQC keeps Cr(3)/Cr(6) in redox equilibrium at the fixed pe, so the
# Python side runs in coupled-redox mode here (SCEPTER-TE itself decouples them).
CASES.append(('edge_CrVI', base(o2_si=-0.68, solution={'Na': 0.01, 'Cl': 0.011, 'Cr(6)': 1e-6}),
              {'metals': ['Cr'], 'redox': 'coupled', 'free': 'CrO4-2'}))
CASES.append(('edge_CrIII', base(pe=2, solution={'Na': 0.01, 'Cl': 0.011, 'Cr(3)': 1e-6}),
              {'metals': ['Cr'], 'redox': 'coupled', 'free': 'Cr(OH)2+'}))
# Arsenic: As(V) in decoupled mode at oxic pe; As(III) coupled at low pe (as Cr(III))
CASES.append(('edge_AsV', base(solution={'Na': 0.01, 'Cl': 0.011, 'As(5)': 1e-6}),
              {'metals': ['As(5)'], 'free': 'H2AsO4-'}))
CASES.append(('edge_AsIII', base(pe=0, solution={'Na': 0.01, 'Cl': 0.011, 'As(3)': 1e-6}),
              {'metals': ['As'], 'redox': 'coupled', 'free': 'H3AsO3'}))
# Mn oxide (HMO) edges for the cations it binds, alone and with Hfo
for m in ['Ni', 'Co', 'Cu', 'Zn', 'Cd', 'Pb']:
    # Co at pe 6: at pe 12 PHREEQC oxidises some Co to Co(III), which SCEPTER-TE does not carry
    CASES.append(('hmo_' + m, base(hfo=None, hmo=HMO, pe=6 if m == 'Co' else 12,
                                   solution={'Na': 0.01, 'Cl': 0.011, m: 1e-6}),
                  {'metals': [m]}))
# Humic acid (simplified WHAM VII, non-electrostatic): 0.1 g HA, alone
for m in ['Ni', 'Co', 'Cu', 'Zn', 'Cd', 'Pb']:
    CASES.append(('hum_' + m, base(hfo=None, hum=0.1, pe=6 if m == 'Co' else 12,
                                   solution={'Na': 0.01, 'Cl': 0.011, m: 1e-6}),
                  {'metals': [m]}))
CASES.append(('hum_soil', base(hfo=None, hum=0.1, pco2_si=-2.0, cec=0.05,
                               solution={'Na': 1e-3, 'K': 2e-4, 'Ca': 2e-3, 'Mg': 5e-4, 'Cl': 4e-3,
                                         'Al': 1e-6, 'Ni': 1e-6, 'Cu': 1e-6, 'Zn': 1e-6,
                                         'Cd': 1e-7, 'Pb': 1e-7},
                               pH=[4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]),
              {'metals': METALS + ['Al']}))
# As(V) with competing phosphate, silicate (basalt) and bicarbonate (pCO2 1%)
CASES.append(('AsV_competition', base(pco2_si=-2.0, solution={
    'Na': 0.01, 'Cl': 0.011, 'Ca': 1e-3, 'As(5)': 1e-6, 'P': 2e-5, 'Si': 5e-4}),
    {'metals': ['As(5)', 'P']}))
CASES.append(('edge_Ni_W0.5', base(W=0.5, solution={'Na': 0.01, 'Cl': 0.011, 'Ni': 1e-6}),
              {'metals': ['Ni']}))
soil = {'Na': 1e-3, 'K': 2e-4, 'Ca': 2e-3, 'Mg': 5e-4, 'Cl': 4e-3, 'S(6)': 5e-4,
        'Ni': 1e-6, 'Cu': 1e-6, 'Zn': 1e-6, 'Cd': 1e-7, 'Pb': 1e-7}
CASES.append(('soil_25C', base(pco2_si=-2.0, cec=0.05, solution=soil),
              {'metals': METALS}))
CASES.append(('soil_10C', base(T=10, pco2_si=-2.0, cec=0.05, solution=soil),
              {'metals': METALS}))
CASES.append(('soil_hfo_hmo', base(pco2_si=-2.0, cec=0.05, hmo=HMO,
                                   solution=dict(soil, Co=1e-7, Mn=1e-5)),
              {'metals': METALS + ['Co', 'Mn']}))
# acid soil: Al on the exchanger, gibbsite cap (Sabah, Vermont)
acid = {'Na': 2e-4, 'K': 1e-4, 'Ca': 5e-4, 'Mg': 1.25e-4, 'Cl': 2e-3, 'Al': 2e-4,
        'Ni': 1e-6, 'Cu': 1e-6}
CASES.append(('acid_soil_Al', base(pco2_si=-2.0, cec=0.05, solution=acid,
                                   minerals=['Gibbsite'], pH=[4.0, 4.3, 4.6, 4.9, 5.2, 5.5, 6.0]),
              {'metals': ['Ni', 'Cu', 'Al']}))
caps = {'Na': 1e-3, 'Cl': 1e-3, 'Cu': 1e-4, 'Zn': 1e-4, 'Pb': 1e-4}
CASES.append(('caps', base(hfo=None, pco2_si=-2.0, solution=caps,
                           minerals=['Tenorite', 'Malachite', 'Cerussite',
                                     'Hydrocerussite', 'Zincite', 'Smithsonite']),
              {'metals': ['Cu', 'Zn', 'Pb']}))


def run_case(name, case, opt):
    labels = list(case['solution'])
    els = sorted({thermo_base(l) for l in labels} | {'Na', 'Cl'}
                 | ({'C'} if case['pco2_si'] is not None else set()))
    val = {'C': ['C(4)'], 'S': ['S(6)']}
    for m in METALS:
        if m == 'Cu':
            val['Cu'] = ['Cu(2)']
    if opt.get('redox', 'decoupled') == 'decoupled':
        val.update({'Co': ['Co(2)'], 'Mn': ['Mn(2)']})
        for l in labels:
            if '(' in l and thermo_base(l) not in val:
                val[thermo_base(l)] = [l]
    sy = system.build(db, els, redox=opt.get('redox', 'decoupled'), valences=val,
                      exchangers=['X'] if case['cec'] else [],
                      surfaces=(['Hfo_w', 'Hfo_s'] if case['hfo'] else [])
                      + (['Hmo_x', 'Hmo_y'] if case.get('hmo') else [])
                      + (list(humic.PER_G) if case.get('hum') else []),
                      phases=case['minerals'] + (['CO2(g)'] if case['pco2_si'] is not None else []))
    species = sy.aq.names + sy.ex.names + sy.surf.names
    rows = pr.run(dict(case, elements=els), species)
    sv = solver.Speciator(sy)
    out = []
    for r in rows:
        tot = {}
        for e in els:
            tot[e] = r['SYS_' + e]
        kw = {}
        if case['hfo']:
            h = case['hfo']
            kw['surfaces'] = {'Hfo': solver.Surface({'Hfo_w': h['weak'], 'Hfo_s': h['strong']},
                                                    h['area'], h['mass'])}
        if case.get('hmo'):
            h = case['hmo']
            kw.setdefault('surfaces', {})['Hmo'] = solver.Surface(
                {'Hmo_x': h['x'], 'Hmo_y': h['y']}, h['area'], h['mass'])
        if case.get('hum'):
            kw.setdefault('surfaces', {})['Hum'] = humic.surface(case['hum'])
        if case['cec']:
            kw['cec'] = {'X': case['cec']}
        if case['pco2_si'] is not None:
            kw['fix_si'] = [('CO2(g)', case['pco2_si'], 'CO3-2')]
            tot.pop('C')
        res = sv.solve(tot, pH=r['pH'], T=case['T'] + 273.15, W=r['W'],
                       pe=r['pe'] if opt.get('redox') == 'coupled' else None,
                       minerals=case['minerals'], **kw)
        for m in opt['metals']:
            lab = m
            k = sy.index(sy.aliases[lab])
            total = r['SYS_' + thermo_base(m)]
            # PHREEQC pools from its own species, using the shared stoichiometry
            # (MOL() returns mol/kgw for exchange and surface species too)
            ph = {'dissolved': sum(r['MOL_' + n] * r['W'] * sy.aq.nu[i, k]
                                   for i, n in enumerate(sy.aq.names)),
                  'exchange': sum(r['MOL_' + n] * r['W'] * sy.ex.nu[i, k]
                                  for i, n in enumerate(sy.ex.names)),
                  'sorbed': sum(r['MOL_' + n] * r['W'] * sy.surf.nu[i, k]
                                for i, n in enumerate(sy.surf.names)),
                  'mineral': sum(r['EQUI_' + p] * sy.phases[p][0][k] for p in case['minerals'])}
            free_sp = opt.get('free', sy.basis[k])
            ph_free = r['MOL_' + free_sp] * r['W']
            py = res.partition[lab]
            row = dict(case=name, T=case['T'], pH=round(r['pH'], 3), metal=lab,
                       total=total, mu_phreeqc=r['mu'], mu_py=res.mu,
                       free_phreeqc=ph_free, free_py=res.aq[free_sp],
                       closure_phreeqc=sum(ph.values()) / total - 1)
            for p in ('dissolved', 'exchange', 'sorbed', 'mineral'):
                row['f_%s_phreeqc' % p] = ph[p] / total
                row['f_%s_py' % p] = py[p] / total
            out.append(row)
    return out


def thermo_base(label):
    return system._base_elem(label)


def main():
    os.makedirs(OUT, exist_ok=True)
    allrows = []
    for name, case, opt in CASES:
        print('running', name, flush=True)
        allrows += run_case(name, case, opt)
    keys = list(allrows[0].keys())
    with open(os.path.join(OUT, 'validation.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(allrows)
    # summary
    lines = ['| case | metal | max abs error in pool fraction | max rel. error free ion | '
             'max rel. error ionic strength | PHREEQC mass closure |',
             '|---|---|---|---|---|---|']
    worst = 0.0
    for name, _, opt in CASES:
        for m in sorted({r['metal'] for r in allrows if r['case'] == name}):
            rr = [r for r in allrows if r['case'] == name and r['metal'] == m]
            df = max(abs(r['f_%s_phreeqc' % p] - r['f_%s_py' % p]) for r in rr
                     for p in ('dissolved', 'exchange', 'sorbed', 'mineral'))
            dfree = max(abs(r['free_py'] / r['free_phreeqc'] - 1) for r in rr)
            dmu = max(abs(r['mu_py'] / r['mu_phreeqc'] - 1) for r in rr)
            clo = max(abs(r['closure_phreeqc']) for r in rr)
            worst = max(worst, df)
            lines.append(f'| {name} | {m} | {df:.1e} | {dfree:.1e} | {dmu:.1e} | {clo:.1e} |')
    md = '\n'.join(lines)
    with open(os.path.join(OUT, 'validation.md'), 'w') as f:
        f.write(md + '\n')
    print(md)
    plot(allrows)
    return worst


def plot(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    panels = [('edge_Ni', 'Ni'), ('edge_Cu', 'Cu'), ('edge_Zn', 'Zn'), ('edge_Cd', 'Cd'),
              ('edge_Pb', 'Pb'), ('edge_CrVI', 'Cr'), ('edge_CrIII', 'Cr'),
              ('soil_10C', 'Ni'), ('soil_10C', 'Cu'), ('caps', 'Cu'), ('caps', 'Zn'), ('caps', 'Pb')]
    fig, axs = plt.subplots(3, 4, figsize=(13, 9), sharex=True, sharey=True)
    cols = {'dissolved': '#2a6fdb', 'sorbed': '#d9822b', 'exchange': '#3a9d5d', 'mineral': '#8e44ad'}
    for ax, (case, m) in zip(axs.flat, panels):
        rr = [r for r in rows if r['case'] == case and r['metal'] == m]
        ph = [r['pH'] for r in rr]
        for p, c in cols.items():
            yp = [r['f_%s_py' % p] for r in rr]
            if max(yp) < 1e-3 and max(r['f_%s_phreeqc' % p] for r in rr) < 1e-3:
                continue
            ax.plot(ph, yp, '-', color=c, label=p + ' (Python)')
            ax.plot(ph, [r['f_%s_phreeqc' % p] for r in rr], 'o', mfc='none', color=c,
                    label=p + ' (PHREEQC)')
        ax.set_title(f'{case}: {m}', fontsize=9)
        ax.set_ylim(-0.03, 1.03)
    for ax in axs[-1]:
        ax.set_xlabel('pH')
    for ax in axs[:, 0]:
        ax.set_ylabel('fraction of total')
    h, l = [], []
    for ax in axs.flat:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh)
                l.append(ll)
    fig.legend(h, l, loc='lower center', ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(OUT, 'pH_edges.png'), dpi=120)


if __name__ == '__main__':
    w = main()
    print(f'worst pool-fraction error: {w:.2e}')
    sys.exit(0 if w < 0.01 else 1)
