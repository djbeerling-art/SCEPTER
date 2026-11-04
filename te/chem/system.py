"""Build a reduced chemical system (components x species) from a Database.

Every species is rewritten in terms of a basis: primary master species
(Ni+2, CO3-2, ...), H+, e-, H2O, exchange masters (X-) and surface masters
(Hfo_wOH, Hfo_sOH). With ``redox='coupled'`` all valence states are linked
through e- (the solver then needs a fixed pe). With ``redox='decoupled'``
the listed valence states, e.g. Cr(3) and Cr(6), become separate components
and no species carries e-; this is the mode SCEPTER-TE uses, since Cr redox
is kinetic (Mn oxides, Fe(II), SOM) rather than equilibrium.
"""
import re
from dataclasses import dataclass, field

import numpy as np

from . import thermo

FIXED = ('H+', 'e-', 'H2O')


def _ex_gamma(db, r, exm):
    """PHREEQC's activity coefficient for an exchange species with -gamma
    (model.cpp calc_gammas): log g = coef*(-A z^2 sqrt(mu)/(1 + a B sqrt(mu))
    + b mu), or Davies (-coef z^2 A (sqrt(mu)/(1+sqrt(mu)) - 0.3 mu)) when
    a = b = 0, where coef and z belong to the LAST aqueous token (not H2O) of
    the reaction rewritten in master species. For CaX2 that is Ca+2 (coef 1,
    z 2); for AlOHX2 it is H+ (coef -1), so its log g is positive.
    Returns (a, b, coef, z, davies)."""
    primaries = set(db.masters.values()) | {'H+', 'H2O', 'e-'}

    def expand(sp, c, depth=0):
        if sp in primaries or sp not in db.species or depth > 5:
            return [(sp, c)]
        rr = db.species[sp]
        out = []
        for t, v in rr.lhs:
            out += expand(t, c * v, depth + 1)
        for t, v in rr.rhs:
            if t != sp:
                out += expand(t, -c * v, depth + 1)
        return out

    toks = []
    for sp, c in r.lhs:
        if sp not in exm:
            toks += expand(sp, c)
    for sp, c in r.rhs:
        if sp != r.defined and sp not in exm:
            toks += expand(sp, -c)
    toks = [(t, c) for t, c in toks if t not in ('H2O', 'e-')]
    coef, zt = (toks[-1][1], thermo.charge(toks[-1][0])) if toks else (0.0, 0.0)
    a, b = r.gamma
    davies = a == 0 and b == 0
    return (a, b, coef, zt, davies)


def _base_elem(e):
    return re.sub(r'\(.*\)$', '', e)


@dataclass
class Group:
    """Species of one kind (aq, ex, surf) with stoichiometry on the basis."""
    names: list
    nu: np.ndarray                 # (n_species, n_basis)
    rxn: list                      # (logk25-shift-ready) list of (Reaction chain, sign)
    z: np.ndarray                  # charge (aq, surf) or number of X (ex)
    gamma: list                    # (a, b) or None
    logk_terms: list               # list of [(Reaction, factor)] for log K(T)

    def logk(self, T):
        return np.array([sum(f * r.logk(T) for r, f in terms) for terms in self.logk_terms])


@dataclass
class System:
    basis: list                    # basis species names
    comp_elem: dict                # basis species -> element/valence label
    aq: Group
    ex: Group
    surf: Group
    phases: dict = field(default_factory=dict)   # name -> (nu vector, logk_terms)
    aliases: dict = field(default_factory=dict)  # 'Cr(6)' -> 'CrO4-2', ...

    def index(self, name):
        return self.basis.index(name)

    def phase_logk(self, name, T):
        return sum(f * r.logk(T) for r, f in self.phases[name][1])


class _Reducer:
    def __init__(self, db):
        self.db = db
        self.cache = {}
        allsp = {}
        allsp.update(db.species)
        allsp.update(db.ex_species)
        allsp.update(db.surf_species)
        self.allsp = allsp

    def is_basis(self, name):
        r = self.allsp.get(name)
        return r is not None and len(r.lhs) == 1 and len(r.rhs) == 1 \
            and r.lhs[0][0] == name and r.rhs[0][0] == name

    def expr(self, name):
        """log a(name) = sum(f*logK) + nu . log a(basis); returns (nu dict, terms)."""
        if name in self.cache:
            return self.cache[name]
        if self.is_basis(name):
            out = ({name: 1.0}, [])
        else:
            r = self.allsp[name]
            cD = dict(r.rhs)[name]
            nu, terms = {}, [(r, 1.0 / cD)]
            for sp, c in r.lhs:
                self._add(nu, terms, sp, c / cD)
            for sp, c in r.rhs:
                if sp != name:
                    self._add(nu, terms, sp, -c / cD)
            out = (nu, terms)
        self.cache[name] = out
        return out

    def _add(self, nu, terms, sp, c):
        n2, t2 = self.expr(sp)
        for k, v in n2.items():
            nu[k] = nu.get(k, 0.0) + c * v
        terms.extend((r, f * c) for r, f in t2)

    def phase_expr(self, phase):
        """SI = nu . log a - logK; logK expressed via terms (sign folded in)."""
        r = self.db.phases[phase]
        nu, terms = {}, [(r, -1.0)]
        # log IAP = sum rhs - sum lhs (excluding the phase formula itself)
        for sp, c in r.rhs:
            self._add(nu, terms, sp, c)
        for sp, c in r.lhs[1:]:
            self._add(nu, terms, sp, -c)
        # terms now carry +logK of the species formation; SI = nu.la + sum(terms)
        return nu, terms


def build(db, elements, redox='decoupled', valences=None, exchangers=(),
          surfaces=(), phases=(), drop_species=()):
    """Assemble a System.

    elements : element names ('Na', 'Ni', 'C', 'Cr', ...).
    valences : for redox='decoupled', the valence states to carry per redox
               element, e.g. {'Cr': ['Cr(3)', 'Cr(6)'], 'C': ['C(4)']}. An
               element not listed is carried in the state of its primary master.
    exchangers / surfaces : exchange master names ('X') / surface masters
               ('Hfo_w', 'Hfo_s') to include.
    phases : phases whose SI the solver may need (gases, minerals).
    """
    red = _Reducer(db)
    valences = valences or {}
    prim = {}   # primary master species -> element
    for el, sp in db.masters.items():
        if '(' not in el and el in elements:
            prim[sp] = el
    # secondary masters for decoupled valence states
    sec = {}    # label -> (master species, nu over basis, terms)
    if redox == 'decoupled':
        for el, labels in valences.items():
            if el not in elements:
                continue
            for lab in labels:
                sp = db.masters[lab]
                if sp in prim:
                    continue
                nu, terms = red.expr(sp)
                sec[lab] = (sp, nu, terms)
    basis = list(prim) + [v[0] for v in sec.values()] + list(FIXED)
    comp_elem = {sp: el for sp, el in prim.items()}
    for lab, (sp, _, _) in sec.items():
        comp_elem[sp] = lab
    exm = [db.ex_masters[x] for x in exchangers]
    sfm = [db.surf_masters[s] for s in surfaces]
    basis += exm + sfm
    for x, sp in zip(exchangers, exm):
        comp_elem[sp] = x
    for s, sp in zip(surfaces, sfm):
        comp_elem[sp] = s
    allowed_full = set(prim) | set(FIXED) | set(exm) | set(sfm)
    # which primary is the valence-reference of each secondary
    sec_by_prim = {}
    for lab, (sp, nu, terms) in sec.items():
        p = [k for k in nu if k in prim]
        assert len(p) == 1, lab
        sec_by_prim.setdefault(p[0], []).append((lab, sp, nu, terms))
    keep_prim_state = {p: True for p in prim}
    for el, labels in valences.items():
        for p, e in prim.items():
            if e == el and redox == 'decoupled':
                keep_prim_state[p] = db.masters[el] in [db.masters[l] for l in labels]

    def rebase(nu, terms):
        """Map a fully reduced expression onto the decoupled basis; None = drop."""
        if any(k not in allowed_full for k in nu if abs(nu[k]) > 1e-12):
            return None
        nu = dict(nu)
        terms = list(terms)
        if redox == 'coupled':
            return nu, terms
        for p in prim:
            n = nu.get(p, 0.0)
            if abs(n) < 1e-12:
                continue
            ne = nu.get('e-', 0.0)
            if abs(ne) < 1e-12 or not sec_by_prim.get(p):
                # no electrons left, or this element carries no secondary
                # valence (the e- belongs to another element's master)
                if not keep_prim_state[p]:
                    return None
                continue
            hit = None
            for lab, sp, snu, sterms in sec_by_prim.get(p, []):
                if abs(snu.get('e-', 0.0) * n - ne) < 1e-9:
                    hit = (sp, snu, sterms)
            if hit is None:
                return None
            sp, snu, sterms = hit
            for k, v in snu.items():
                nu[k] = nu.get(k, 0.0) - n * v
            nu[sp] = nu.get(sp, 0.0) + n
            terms += [(r, -n * f) for r, f in sterms]
        if abs(nu.get('e-', 0.0)) > 1e-9:
            return None
        nu.pop('e-', None)
        return nu, terms

    def group(names, kind):
        rows, kept, terms_l, z, gam = [], [], [], [], []
        for name in names:
            if name in drop_species:
                continue
            res = rebase(*red.expr(name))
            if res is None:
                continue
            nu, terms = res
            row = np.zeros(len(basis))
            for k, v in nu.items():
                if abs(v) > 1e-12:
                    row[basis.index(k)] = v
            if kind == 'aq' and name in ('H2O',):
                continue
            rows.append(row)
            kept.append(name)
            terms_l.append(terms)
            r = red.allsp[name]
            if kind == 'ex':
                z.append(sum(row[basis.index(m)] for m in exm))
            else:
                z.append(thermo.charge(name))
            if kind == 'ex' and r.gamma is not None:
                gam.append(_ex_gamma(db, r, exm))
            else:
                gam.append(r.gamma)
        nu = np.array(rows) if rows else np.zeros((0, len(basis)))
        return Group(kept, nu, None, np.array(z), gam, terms_l)

    aq_names = [s for s in db.species if s not in ('e-',)]
    aq = group(aq_names, 'aq')
    ex = group([s for s in db.ex_species if s not in exm], 'ex')
    surf = group([s for s in db.surf_species if s not in sfm], 'surf')
    # basis exchange / surface masters must themselves be species
    # the free exchanger X- is a master only (no moles, as in PHREEQC)
    for m in sfm:
        _append(surf, m, basis, 0.0)
    ph = {}
    for p in phases:
        res = rebase(*red.phase_expr(p))
        if res is None:
            raise ValueError(f'phase {p} not representable in this basis')
        nu, terms = res
        row = np.zeros(len(basis))
        for k, v in nu.items():
            row[basis.index(k)] = v
        ph[p] = (row, terms)
    aliases = {lab: sp for lab, sp in db.masters.items()
               if sp in basis and _base_elem(lab) in elements}
    # an element whose primary valence is not carried (e.g. minteq's Co+3 when
    # only Co(2) is kept) is labelled by its first carried valence
    for p, el in prim.items():
        if not keep_prim_state[p] and valences.get(el):
            aliases[el] = db.masters[valences[el][0]]
    return System(basis, comp_elem, aq, ex, surf, ph, aliases)


def _append(g, name, basis, z):
    row = np.zeros(len(basis))
    row[basis.index(name)] = 1.0
    g.names.append(name)
    g.nu = np.vstack([g.nu, row])
    g.z = np.append(g.z, z)
    g.gamma.append(None)
    g.logk_terms.append([])
