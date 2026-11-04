"""Humic acid site amounts for the simplified WHAM Model VII surface in
te_extra.dat (Hum_aa-ad type A, Hum_ba-bd type B, Hum_sa/sb strong;
non-electrostatic)."""
NA = 3.4e-3          # mol type-A sites per g HA (Tipping et al. 2011)

PER_G = {**{f'Hum_a{c}': 0.901 * NA / 4 for c in 'abcd'},
         **{f'Hum_b{c}': NA / 8 for c in 'abcd'},
         'Hum_sa': 0.09 * NA, 'Hum_sb': 0.009 * NA}


def sites(g):
    """Surface master label -> mol sites for g grams of humic acid."""
    return {k: v * g for k, v in PER_G.items()}


def surface(g):
    from .solver import Surface
    return Surface(sites(g), area=0.0, mass=g, edl=False)


# Dissolved organic matter (DOM): a mobile copy of the humic surface ('Dom_*'), same
# binding constants. Metal bound to it travels with the porewater and is counted as
# dissolved. Its amount follows DOC, which rises with pH (more humic matter dissolves
# as its acid groups deprotonate).
DOM_PER_G = {k.replace('Hum_', 'Dom_'): v for k, v in PER_G.items()}
OM_PER_DOC = 2.0        # g organic matter per g C


def add_dom(db):
    """Clone the Hum_* surface masters and species in a thermo.Database as Dom_*."""
    import copy
    for k, v in list(db.surf_masters.items()):
        if k.startswith('Hum_'):
            db.surf_masters[k.replace('Hum_', 'Dom_')] = v.replace('Hum_', 'Dom_')
    for k, rx in list(db.surf_species.items()):
        if k.startswith('Hum_'):
            r = copy.deepcopy(rx)
            r.lhs = [(n.replace('Hum_', 'Dom_'), c) for n, c in r.lhs]
            r.rhs = [(n.replace('Hum_', 'Dom_'), c) for n, c in r.rhs]
            r.defined = r.defined.replace('Hum_', 'Dom_')
            db.surf_species[k.replace('Hum_', 'Dom_')] = r
    return db


def doc_mg_per_L(pH, doc_ref, pH_ref, slope):
    """DOC (mg C/L) = doc_ref * 10**(slope * (pH - pH_ref))."""
    return doc_ref * 10 ** (slope * (pH - pH_ref))


def dom_surface(g):
    from .solver import Surface
    return Surface({k: v * g for k, v in DOM_PER_G.items()}, area=0.0, mass=g, edl=False)
