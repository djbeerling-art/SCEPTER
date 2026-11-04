"""Clay edge sites: 2-site protolysis non-electrostatic surface complexation
(2SPNE SC) of Bradbury & Baeyens for 2:1 clays (montmorillonite, illite).

Site capacities (mol/kg clay): strong 2.0e-3, weak (W1) 4.0e-2, both with
protolysis log K +4.5 / -7.9 (Bradbury & Baeyens 1997, 2005, 2009).
Metal constants for the first complex, S-OH + M+2 = S-OM+ + H+:
Ni from the Ni data (strong -0.6, weak -3.3; B&B 1997); the others from the
B&B (2005) montmorillonite linear free-energy relation
log SK = 8.1 + 0.90 log OHK1, with first-hydrolysis constants (Baes & Mesmer),
and the weak site 2.7 log units below the strong site, as for Ni. The
hydrolysed surface complexes matter only above pH ~8 and are left out.
"""
from .thermo import Reaction

SITES = {'Clay_s': 2.0e-3, 'Clay_w': 4.0e-2}        # mol per g x 1e3 -> per kg; see sites()
LOG_OHK1 = {'Co': -9.65, 'Zn': -8.96, 'Cd': -10.08, 'Pb': -7.71, 'Cu': -7.95}
LOGK_S = {'Ni': -0.6, **{m: round(8.1 + 0.90 * v, 2) for m, v in LOG_OHK1.items()}}
DW = 2.7                                             # strong minus weak


# Measured 2SPNE constants for montmorillonite (strong, weak W1), ClaySor 2023 database
# (Marinich et al. 2025, Appl. Geochem. 191:106510); Zn weak from the same table.
CLAYSOR_MONT = {'Zn': (0.98, -2.36), 'Pb': (1.37, None), 'Co': (-0.18, None), 'Cd': (-1.03, None), 'Ni': (-0.81, None)}


def add_clay(db, overrides=None):
    """Add Clay_s / Clay_w masters and species to a thermo.Database.
    overrides: {metal: (log K strong, log K weak or None -> strong - DW)} replacing the defaults."""
    ks = {m: (k, k - DW) for m, k in LOGK_S.items()}
    for m, (kst, kw) in (overrides or {}).items():
        ks[m] = (kst, kst - DW if kw is None else kw)
    for s in SITES:
        oh = f'{s}OH'
        db.surf_masters[s] = oh
        db.surf_species[oh] = Reaction([(oh, 1.0)], [(oh, 1.0)], 0.0, defined=oh)
        db.surf_species[f'{s}OH2+'] = Reaction([(oh, 1.0), ('H+', 1.0)], [(f'{s}OH2+', 1.0)], 4.5, defined=f'{s}OH2+')
        db.surf_species[f'{s}O-'] = Reaction([(oh, 1.0)], [(f'{s}O-', 1.0), ('H+', 1.0)], -7.9, defined=f'{s}O-')
        for m, (kst, kw) in ks.items():
            sp = f'{s}O{m}+'
            db.surf_species[sp] = Reaction([(oh, 1.0), (f'{m}+2', 1.0)], [(sp, 1.0), ('H+', 1.0)],
                                           kw if s == 'Clay_w' else kst, defined=sp)
    return db


def surface(g):
    """2SPNE surface for g grams of 2:1 clay (non-electrostatic)."""
    from .solver import Surface
    return Surface({k: v * g / 1000.0 for k, v in SITES.items()}, area=0.0, mass=g, edl=False)
