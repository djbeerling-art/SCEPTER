"""Virtual soil extractions: re-equilibrate a modelled soil cell with a lab
extractant, so that model and measurement report the same operational pool.

The cell's metal and host totals (exchanger, oxide and humic sites, porewater)
are put into the extractant volume at the protocol's solution:soil ratio, the
pH is fixed at the extractant pH, CO2 degasses to air, and whatever is
dissolved at equilibrium is the modelled extractable metal (mg/kg soil).

Caveat: activity corrections use the solver's Davies/extended Debye-Hueckel
terms, which are rough at ~1 M ionic strength; ratios (basalt/control) are
much less sensitive to this than absolute values.
"""
import numpy as np

# reagent totals in mol/L of extractant; ratio = L extractant per kg soil
PROTOCOLS = {
    # 1 M ammonium acetate, pH 7 (Sabah, Corn Belt)
    'NH4OAc': dict(pH=7.0, ratio=10.0, reagents={'N(-3)': 1.0, 'Acetate': 1.0}),
    # Modified Morgan: 0.62 M NH4OH + 1.25 M acetic acid, pH 4.8 (Vermont)
    'MMorgan': dict(pH=4.8, ratio=5.0, reagents={'N(-3)': 0.62, 'Acetate': 1.25}),
    # Mehlich-3: 0.2 M HOAc, 0.25 M NH4NO3, 0.015 M NH4F, 0.013 M HNO3, 0.001 M EDTA, pH 2.5
    'Mehlich3': dict(pH=2.5, ratio=10.0, reagents={'Acetate': 0.2, 'N(-3)': 0.265,
                                                    'N(5)': 0.263, 'F': 0.015, 'Edta': 0.001}),
}
ELEMENTS = {'NH4OAc': ['N', 'Acetate'], 'MMorgan': ['N', 'Acetate'],
            'Mehlich3': ['N', 'Acetate', 'F', 'Edta']}
VALENCES = {'N': ['N(-3)', 'N(5)']}


def extract(sv, cell_totals, *, soil_kg, cec, surfaces, protocol, T=293.15, guess=None):
    """cell_totals: {label: mol in the cell} (host majors + metals, as in the
    TE module); soil_kg: soil mass of the cell; cec: exchanger charge (mol);
    surfaces: {name: Surface} for the cell. Returns (result, W)."""
    p = PROTOCOLS[protocol]
    W = p['ratio'] * soil_kg                     # L of extractant
    tot = dict(cell_totals)
    for lab, c in p['reagents'].items():
        tot[lab] = tot.get(lab, 0.0) + c * W
    tot.pop('C', None)
    r = sv.solve(tot, pH=p['pH'], W=W, T=T, cec={'X': cec}, surfaces=surfaces,
                 fix_si=[('CO2(g)', -3.4, 'CO3-2')], guess=guess)
    return r, W


def extracted_cells(sv, res, snap, k, cells, metals, *, protocol, soil_kg_per_m3):
    """Per-cell extracted metal (mol per m2 of cell) and soil mass (kg per m2)
    at snapshot k of a scepter_te.TEResult (res.host, res.surfaces, res.cec_x, res.pools)."""
    out = np.zeros((len(cells), len(metals))); soil = np.zeros(len(cells))
    cec = res.cec_x(snap)
    for n, i in enumerate(cells):
        vol = snap.dz[i]                                     # m3 per m2
        sk = soil_kg_per_m3[i] * vol
        tot = {l: v * vol for l, v in res.host[k][i].items()}
        for j, m in enumerate(metals):
            tot[res.labels[m]] = max(res.pools['total'][k][i, j], 1e-30) * vol
        surf = {nm: _scale(s, vol) for nm, s in res.surfaces(i).items()}
        r, W = extract(sv, tot, soil_kg=sk, cec=float(cec[i]) * vol, surfaces=surf,
                       protocol=protocol)
        for j, m in enumerate(metals):
            out[n, j] = r.partition[res.labels[m]]['dissolved']
        soil[n] = sk
    return out, soil


def extractable_mgkg(sv, res, snap, k, cells, metals, mw, *, protocol, soil_kg_per_m3):
    """Mass-weighted extractable metal (mg/kg) over `cells` at snapshot k."""
    out, soil = extracted_cells(sv, res, snap, k, cells, metals, protocol=protocol,
                                soil_kg_per_m3=soil_kg_per_m3)
    mwv = np.array([mw[m] for m in metals])
    return (out * mwv * 1e3).sum(0) / soil.sum()


def _scale(s, f):
    """Surface sites for a cell are per m3 bulk in the TE module; scale to the cell."""
    from .solver import Surface
    return Surface({k: v * f for k, v in s.sites.items()}, s.area, s.mass * f, s.edl)
