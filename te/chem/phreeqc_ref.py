"""Run the same batch problems in PHREEQC (IPhreeqc via phreeqpython).

Each case fixes pH (Fix_H+ with NaOH), optionally pe (Fix_pe with O2(g)) and
pCO2 (CO2(g)), and equilibrates a solution with an exchanger, Hfo surface and
optional minerals. PHREEQC's system totals after the acid/base additions are
returned so the Python solver can be driven with exactly the same inputs.
"""
import os

from phreeqpython.viphreeqc import VIPhreeqc

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(HERE, 'data', 'minteq.v4.dat')
EXTRA = os.path.join(HERE, 'data', 'te_extra.dat')


def _ip():
    ip = VIPhreeqc()
    ip.load_database(DB)
    return ip


def _surf_lines(case):
    out = []
    if case.get('hfo'):
        h = case['hfo']
        out += [f"  Hfo_wOH {h['weak']} {h['area']} {h['mass']}", f"  Hfo_sOH {h['strong']}"]
    if case.get('hmo'):
        h = case['hmo']
        out += [f"  Hmo_xOH {h['x']} {h['area']} {h['mass']}", f"  Hmo_yOH {h['y']}"]
    if case.get('hum'):
        from .humic import sites
        out += [f"  {k}H {v}" for k, v in sites(case['hum']).items()]
        out.append('  -no_edl')   # applies to the whole SURFACE block
    return out


def run(case, species):
    """case: dict with T (C), W (kg), solution {label: mol/kgw}, pH list,
    pe (or None; fixed with Fix_pe), o2_si (or None; pe from O2(g) instead),
    pco2_si (or None), cec (mol), hfo (dict weak, strong, mass,
    area) or None, hmo (dict x, y, mass, area) or None, minerals list, elements list (for SYS output).
    species: aqueous + exchange + surface species names to punch.
    Returns a list of dict rows, one per pH."""
    extra = open(EXTRA).read()
    sol = '\n'.join(f'  {k} {v}' for k, v in case['solution'].items())
    lines = [extra, 'SOLUTION 1', f"  temp {case.get('T', 25)}", '  units mol/kgw',
             f"  -water {case.get('W', 1.0)}", '  pH 3.0', sol]
    if case.get('pe') is not None:
        lines.append(f"  pe {case['pe']}")
    lines.append('END')
    names = ['pH', 'pe', 'mu', 'W']
    els = case['elements']
    punch = ['10 PUNCH -LA("H+"), -LA("e-"), MU, TOT("water")']
    ln = 20
    for e in els:
        punch.append(f'{ln} PUNCH SYS("{e}")')
        names.append('SYS_' + e)
        ln += 1
    for sp in species:
        punch.append(f'{ln} PUNCH MOL("{sp}"), LG("{sp}")')
        names += ['MOL_' + sp, 'LG_' + sp]
        ln += 1
    for m in case.get('minerals', []):
        punch.append(f'{ln} PUNCH EQUI("{m}")')
        names.append('EQUI_' + m)
        ln += 1
    head = ' '.join(f'"{n}"' for n in names)
    rows = []
    for pH in case['pH']:
        b = list(lines)
        b += ['USE SOLUTION 1', 'EQUILIBRIUM_PHASES 1',
              f'  Fix_H+ {-pH} NaOH 10']
        if case.get('o2_si') is not None:
            b.append(f"  O2(g) {case['o2_si']} 10")
        elif case.get('pe') is not None:
            b.append(f"  Fix_pe {-case['pe']} O2(g) 10")
        if case.get('pco2_si') is not None:
            b.append(f"  CO2(g) {case['pco2_si']} 10")
        for m in case.get('minerals', []):
            b.append(f'  {m} 0 0')
        if case.get('cec'):
            b += ['EXCHANGE 1', f"  X {case['cec']}", '  -equilibrate 1']
        if case.get('hfo') or case.get('hmo') or case.get('hum'):
            b += ['SURFACE 1'] + _surf_lines(case) + ['  -equilibrate 1']
        b += ['USER_PUNCH', '-headings ' + head] + punch + \
             ['SELECTED_OUTPUT', '-reset false', 'END']
        ip = _ip()
        ip.run_string('\n'.join(b))
        arr = ip.get_selected_output_array()
        row = dict(zip(names, arr[-1]))
        rows.append(row)
    return rows


def run_transport(spec):
    """1-D advective column in PHREEQC TRANSPORT (no dispersion, Courant 1).
    spec: ncell, shifts, inflow {label: mol/kgw, 'pH':, 'charge':},
    soil {same}, cec (mol per cell), hfo dict or None, elements (SYS
    output), metals, ex_species {metal: [exchange species]}.
    Returns (init_row, rows) where rows[(shift, cell)] = dict."""
    extra = open(EXTRA).read()

    def sol(n, s):
        out = [f'SOLUTION {n}', '  units mol/kgw', f"  pH {s['pH']}", '  pe 12']
        for k, v in s.items():
            if k in ('pH',):
                continue
            out.append(f'  {k} {v}' + (' charge' if k == s.get('charge_on') else ''))
        return [l for l in out if not l.strip().startswith('charge_on')]

    n = spec['ncell']
    L = [extra]
    L += sol(0, spec['inflow'])
    L += sol(f'1-{n}', spec['soil'])
    L += [f'EXCHANGE 1-{n}', f"  X {spec['cec']}", '  -equilibrate 1']
    if spec.get('hfo') or spec.get('hmo'):
        L += [f'SURFACE 1-{n}'] + _surf_lines(spec) + ['  -equilibrate 1']
    if spec.get('phases'):
        L.append(f'EQUILIBRIUM_PHASES 1-{n}')
        L += [f'  {p} {si} {amt}' for p, si, amt in spec['phases']]
    L.append('END')
    names = ['sim', 'step', 'cell', 'pH', 'pe', 'W']
    punch = ['10 PUNCH SIM_NO, STEP_NO, CELL_NO, -LA("H+"), -LA("e-"), TOT("water")']
    ln = 20
    for e in spec['elements']:
        punch.append(f'{ln} PUNCH SYS("{e}")'); names.append('SYS_' + e); ln += 1
    for m in spec['metals']:
        exs = ' + '.join(f'MOL("{x}")' for x in spec['ex_species'][m])
        punch.append(f'{ln} PUNCH TOT("{m}"), {exs}')
        names += ['TOT_' + m, 'EX_' + m]
        ln += 1
    head = ' '.join(f'"{x}"' for x in names)
    sel = ['USER_PUNCH', '-headings ' + head] + punch + \
          ['SELECTED_OUTPUT', '-reset false', '-high_precision true']
    # batch re-equilibration of cell 1 gives the initial state (all cells equal)
    L += sel + ['USE solution 1', 'USE exchange 1']
    if spec.get('hfo') or spec.get('hmo'):
        L.append('USE surface 1')
    if spec.get('phases'):
        L.append('USE equilibrium_phases 1')
    L.append('END')
    L += ['TRANSPORT', f'  -cells {n}', f"  -shifts {spec['shifts']}",
          '  -dispersivity 0', '  -diffusion_coefficient 0',
          f'  -punch_cells 1-{n}', '  -punch_frequency 1', 'END']
    ip = _ip()
    ip.run_string('\n'.join(L))
    arr = ip.get_selected_output_array()
    hdr, data = arr[0], arr[1:]
    init, rows = None, {}
    for a in data:
        r = dict(zip(names, a))
        if int(r['sim']) == 2:
            init = r
        elif int(r['sim']) == 3 and int(r['step']) > 0:
            rows[(int(r['step']), int(r['cell']))] = r
    return init, rows
