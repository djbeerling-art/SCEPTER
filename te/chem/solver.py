"""Equilibrium speciation for SCEPTER-TE: aqueous + exchange + DLM surfaces
+ solubility caps, at fixed pH (and pe if redox is coupled).

pH, pe and water mass come from SCEPTER's converged state; every other
component is given as a system total (moles). Optionally a component is set
by a fixed saturation index instead (e.g. CO3-2 by CO2(g) at the soil pCO2).

Conventions follow PHREEQC so the two can be compared species by species:
  aqueous     log m_i = log K_i + sum nu_ik log a_k - log gamma_i
  exchange    beta_i  = K_i prod a_k^nu_ik / gamma_i   (Gaines-Thomas),
              moles_i = beta_i * CEC / z_i
  surface     moles_i = K_i prod a_k^nu_ik exp(-z_i F psi / RT)
              (surface master 'activity' = moles, aqueous = molal)
  DLM         sigma = sqrt(8000 eps eps0 R T I) sinh(F psi / 2RT)
"""
import re
import math
from dataclasses import dataclass, field

import numpy as np

from .system import FIXED

LN10 = math.log(10.0)
F = 96485.33212
R = 8.314462618
EPS0 = 8.8541878128e-12


def water_props(T, P=1.01325):
    """Relative permittivity (Bradley & Pitzer 1979, as PHREEQC) and density
    (Kell 1975) of water at T (K), P (bar)."""
    U = (3.4279e2, -5.0866e-3, 9.4690e-7, -2.0525, 3.1159e3, -1.8289e2,
         -8.0325e3, 4.21452e6, 2.1417)
    e1000 = U[0] * math.exp(U[1] * T + U[2] * T * T)
    C = U[3] + U[4] / (U[5] + T)
    B = U[6] + U[7] / T + U[8] * T
    eps = e1000 + C * math.log((B + P) / (B + 1000.0))
    t = T - 273.15
    rho = (999.83952 + 16.945176 * t - 7.9870401e-3 * t ** 2 - 46.170461e-6 * t ** 3
           + 105.56302e-9 * t ** 4 - 280.54253e-12 * t ** 5) / (1 + 16.879850e-3 * t) / 1000.
    return eps, rho


def dh_ab(T):
    eps, rho = water_props(T)
    A = 1.82483e6 * math.sqrt(rho) / (eps * T) ** 1.5
    B = 50.2916 * math.sqrt(rho) / math.sqrt(eps * T)   # per Angstrom
    return A, B


def log_gamma(group, mu, T, kind):
    A, B = dh_ab(T)
    s = math.sqrt(mu)
    out = np.zeros(len(group.names))
    for i, (z, g) in enumerate(zip(group.z, group.gamma)):
        if kind == 'surf':
            continue
        if kind == 'ex':           # g = (a, b, coef, z, davies), see system._ex_gamma
            if g is None:
                continue
            a, b, coef, zt, davies = g
            if davies:
                out[i] = -coef * zt * zt * A * (s / (1 + s) - 0.3 * mu)
            else:
                out[i] = coef * (-A * zt * zt * s / (1 + B * a * s) + b * mu)
            continue
        if g is not None:        # WATEQ DH (PHREEQC -gamma a b); a = 0 gives
            # the limiting law, and a neutral species gets b*mu
            out[i] = -A * z * z * s / (1 + B * g[0] * s) + g[1] * mu
        elif z == 0:             # PHREEQC default for neutral species
            out[i] = 0.1 * mu
        else:                    # Davies
            out[i] = -A * z * z * (s / (1 + s) - 0.3 * mu)
    return out


@dataclass
class Surface:
    """A surface type sharing one electrostatic plane (e.g. Hfo_w + Hfo_s)."""
    sites: dict                  # surface master label -> moles of sites
    area: float = 600.0          # m2/g
    mass: float = 0.0            # g
    edl: bool = True             # False: non-electrostatic


@dataclass
class Result:
    T: float
    W: float
    pH: float
    mu: float
    la: dict                     # log activity of every basis species
    aq: dict                     # species -> moles
    ex: dict
    surf: dict
    minerals: dict
    psi: dict                    # surface -> potential (V)
    lg: dict                     # species -> log gamma
    iterations: int
    partition: dict = field(default_factory=dict)  # component -> {pool: moles}


class Speciator:
    def __init__(self, system):
        self.s = system

    def solve(self, totals, *, pH, pe=None, T=298.15, W=1.0, cec=None,
              surfaces=None, fix_si=(), minerals=(), tol=1e-11, maxit=200,
              guess=None):
        """totals: {element or valence label: moles}; cec: {'X': moles};
        surfaces: {'Hfo': Surface}; fix_si: [(phase, SI, basis species)];
        minerals: [phase, ...] allowed to precipitate (initially absent)."""
        s = self.s
        nb = len(s.basis)
        cec = cec or {}
        surfaces = surfaces or {}
        lab2sp = {v: k for k, v in s.comp_elem.items()}
        lab2sp.update(s.aliases)
        T_k = np.zeros(nb)
        present = np.zeros(nb, bool)
        for lab, v in list(totals.items()) + list(cec.items()):
            k = s.index(lab2sp[lab])
            T_k[k] = v
            present[k] = v > 0
        surf_of = {}
        for gname, sf in surfaces.items():
            for lab, v in sf.sites.items():
                k = s.index(lab2sp[lab])
                T_k[k] = v
                present[k] = v > 0
                surf_of[k] = gname
        fixed_si = {s.index(sp): (ph, si) for ph, si, sp in fix_si}
        for k in fixed_si:
            present[k] = True
        la = np.full(nb, -300.0)
        la[s.index('H+')] = -pH
        la[s.index('H2O')] = 0.0
        la[s.index('e-')] = -pe if pe is not None else 0.0
        fixed = np.array([b in FIXED for b in s.basis])
        free = np.where(present & ~fixed)[0]
        for k in free:
            la[k] = math.log10(max(T_k[k], 1e-30) / W)
            if s.basis[k] in [lab2sp.get(x) for x in cec]:
                la[k] = 0.0
            if k in surf_of:
                la[k] = math.log10(T_k[k])
        for k, (ph, si) in fixed_si.items():
            la[k] = -5.0
        if guess is not None:
            for k in free:
                if s.basis[k] in guess:
                    la[k] = guess[s.basis[k]]
        gnames = list(surfaces)
        y = np.zeros(len(gnames))
        mins = list(minerals)
        active, nm = [], {}

        # species that can exist (all their basis species present or fixed)
        def mask(g):
            ok = np.all((g.nu == 0) | (present | fixed)[None, :], axis=1)
            return ok
        m_aq, m_ex, m_sf = mask(s.aq), mask(s.ex), mask(s.surf)
        lk_aq, lk_ex, lk_sf = s.aq.logk(T), s.ex.logk(T), s.surf.logk(T)
        exm = {k: lab2sp_k for k, lab2sp_k in
               ((s.index(lab2sp[x]), x) for x in cec)}
        # exchanger of each ex species and site total
        ex_tot = np.zeros(len(s.ex.names))
        for k in exm:
            ex_tot += (s.ex.nu[:, k] != 0) * T_k[k]
        sf_grp = np.full(len(s.surf.names), -1)
        for k, g in surf_of.items():
            sf_grp[s.surf.nu[:, k] != 0] = gnames.index(g)
        edl = np.array([surfaces[g].edl for g in gnames], bool)
        area = np.array([surfaces[g].area * surfaces[g].mass for g in gnames])
        ph_si = {p: (s.phases[p][0], s.phase_logk(p, T)) for p in
                 set(mins) | {v[0] for v in fixed_si.values()}}
        eps, _ = water_props(T)
        mu = 1e-3
        lg_aq = log_gamma(s.aq, mu, T, 'aq')
        lg_ex = log_gamma(s.ex, mu, T, 'ex')

        def species(la, y):
            lm = lk_aq + s.aq.nu @ la - lg_aq
            n_aq = np.where(m_aq, W * 10.0 ** np.minimum(lm, 30), 0.0)
            lb = lk_ex + s.ex.nu @ la - lg_ex
            with np.errstate(divide='ignore', invalid='ignore'):
                n_ex = np.where(m_ex & (ex_tot > 0),
                                10.0 ** np.minimum(lb, 30) * ex_tot / np.where(s.ex.z > 0, s.ex.z, 1), 0.0)
            ls = lk_sf + s.surf.nu @ la
            coul = np.zeros(len(s.surf.names))
            for gi in range(len(gnames)):
                if edl[gi]:
                    coul[sf_grp == gi] = -s.surf.z[sf_grp == gi] * y[gi]
            n_sf = np.where(m_sf & (sf_grp >= 0),
                            10.0 ** np.minimum(ls + coul / LN10, 30), 0.0)
            return n_aq, n_ex, n_sf

        # starting point: fixed-SI components set directly, then a few
        # continued-fraction sweeps on the mass balances (robust far from
        # the solution, e.g. strongly sorbed metals or a fresh exchanger)
        def set_fixed_si():
            for k, (ph, si) in fixed_si.items():
                nu_p, c = ph_si[ph]
                la[k] = 0.0
                la[k] = (si - c - nu_p @ la) / nu_p[k]
        set_fixed_si()
        mb = [k for k in free if k not in fixed_si]
        for _ in range(60):
            n_aq, n_ex, n_sf = species(la, y)
            worst = 0.0
            for k in mb:
                tot = s.aq.nu[:, k] @ n_aq + s.ex.nu[:, k] @ n_ex + s.surf.nu[:, k] @ n_sf
                if tot <= 0:
                    continue
                nmax = max(np.max(np.abs(g.nu[:, k])) for g in (s.aq, s.ex, s.surf) if len(g.names))
                step = math.log10(T_k[k] / tot) / max(nmax, 1.0)
                la[k] += max(min(step, 2.0), -2.0)
                worst = max(worst, abs(step))
            set_fixed_si()
            if worst < 1e-3:
                break

        it_total = 0
        for _outer in range(len(mins) * 2 + 3):
            for it in range(maxit):
                it_total += 1
                n_aq, n_ex, n_sf = species(la, y)
                # update ionic strength, activity coefficients, water activity
                m = n_aq / W
                mu_new = 0.5 * np.sum(m * s.aq.z ** 2)
                mu = max(mu_new, 1e-12)
                lg_aq = log_gamma(s.aq, mu, T, 'aq')
                lg_ex = log_gamma(s.ex, mu, T, 'ex')
                la[s.index('H2O')] = math.log10(max(1 - 0.017 * np.sum(m[s.aq.z != 0.5]), 0.5))
                n_aq, n_ex, n_sf = species(la, y)
                nvar = len(free) + len(gnames) + len(active)
                Fv = np.zeros(nvar)
                J = np.zeros((nvar, nvar))
                scale = np.ones(nvar)
                # mass balances (or fixed-SI rows)
                for r, k in enumerate(free):
                    if k in fixed_si:
                        ph, si = fixed_si[k]
                        nu_p, c = ph_si[ph]
                        Fv[r] = nu_p @ la + c - si
                        J[r, :len(free)] = nu_p[free]
                        continue
                    tot = 0.0
                    for g, n in ((s.aq, n_aq), (s.ex, n_ex), (s.surf, n_sf)):
                        nk = g.nu[:, k]
                        tot += nk @ n
                        J[r, :len(free)] += LN10 * (nk * n) @ g.nu[:, free]
                        if g is s.surf:
                            for gi in range(len(gnames)):
                                if edl[gi]:
                                    sel = sf_grp == gi
                                    J[r, len(free) + gi] += -np.sum(nk[sel] * s.surf.z[sel] * n[sel])
                    for j, p in enumerate(active):
                        tot += ph_si[p][0][k] * nm[p]
                        J[r, len(free) + len(gnames) + j] = ph_si[p][0][k]
                    Fv[r] = tot - T_k[k]
                    scale[r] = max(T_k[k], 1e-30)
                # surface charge balance
                for gi, g in enumerate(gnames):
                    r = len(free) + gi
                    sel = sf_grp == gi
                    z = s.surf.z[sel]
                    n = n_sf[sel]
                    if not edl[gi]:
                        J[r, r] = 1.0
                        continue
                    G = math.sqrt(8000.0 * eps * EPS0 * R * T * mu)
                    Fv[r] = np.sum(z * n) - area[gi] * G / F * math.sinh(y[gi] / 2)
                    J[r, :len(free)] = LN10 * (z * n) @ s.surf.nu[sel][:, free]
                    J[r, r] = -np.sum(z * z * n) - area[gi] * G / F * 0.5 * math.cosh(y[gi] / 2)
                    scale[r] = max(sum(T_k[k] for k, gg in surf_of.items() if gg == g), 1e-30)
                # active minerals: SI = 0
                for j, p in enumerate(active):
                    r = len(free) + len(gnames) + j
                    nu_p, c = ph_si[p]
                    Fv[r] = nu_p @ la + c
                    J[r, :len(free)] = nu_p[free]
                res = np.max(np.abs(Fv / scale)) if nvar else 0.0
                if res < tol and it > 2:
                    break
                Js = J / scale[:, None]
                try:
                    d = np.linalg.solve(Js, -Fv / scale)
                except np.linalg.LinAlgError:
                    d = np.linalg.lstsq(Js, -Fv / scale, rcond=None)[0]
                dl = d[:len(free)]
                dy = d[len(free):len(free) + len(gnames)]
                f = 1.0
                big = np.max(np.abs(dl)) if len(dl) else 0.0
                if big > 1.0:
                    f = 1.0 / big
                if len(dy) and np.max(np.abs(dy)) * f > 2.0:
                    f = 2.0 / np.max(np.abs(dy))
                la[free] += f * dl
                y += f * dy
                for j, p in enumerate(active):
                    nm[p] = nm[p] + f * d[len(free) + len(gnames) + j]
            else:
                raise RuntimeError(f'speciation did not converge (residual {res:.2e})')
            # active-set update for minerals
            changed = False
            for p in list(active):
                if nm[p] < 0:
                    active.remove(p)
                    nm.pop(p)
                    changed = True
            si = {p: ph_si[p][0] @ la + ph_si[p][1] for p in mins if p not in active}
            if si:
                p, v = max(si.items(), key=lambda kv: kv[1])
                if v > 1e-8:
                    # phase rule: if p is not independent of the active set on
                    # the mass-balance components, it replaces the phases it
                    # conflicts with (it is the more supersaturated one)
                    def rank(ps):
                        return np.linalg.matrix_rank(np.array([ph_si[q][0][mb] for q in ps])) if ps else 0
                    if rank(active + [p]) < len(active) + 1:
                        for q in list(active):
                            rest = [a for a in active if a != q]
                            if rank(rest + [p]) == len(rest) + 1:
                                active.remove(q)
                                nm.pop(q)
                                break
                    active.append(p)
                    nm[p] = 0.0
                    changed = True
            if not changed:
                break

        lg = dict(zip(s.aq.names, lg_aq))
        lg.update(zip(s.ex.names, lg_ex))
        res = Result(T=T, W=W, pH=pH, mu=mu,
                     la=dict(zip(s.basis, la)),
                     aq={n: v for n, v, ok in zip(s.aq.names, n_aq, m_aq) if ok},
                     ex={n: v for n, v, ok in zip(s.ex.names, n_ex, m_ex) if ok and v > 0},
                     surf={n: v for n, v, ok in zip(s.surf.names, n_sf, m_sf) if ok and v > 0},
                     minerals={p: nm.get(p, 0.0) for p in mins},
                     psi={g: y[i] * R * T / F for i, g in enumerate(gnames)},
                     lg=lg, iterations=it_total)
        hum_rows = np.array([n.startswith('Hum_') for n in s.surf.names], dtype=bool)   # humic (organic) surface species
        dom_rows = np.array([n.startswith('Dom_') for n in s.surf.names], dtype=bool)   # dissolved organic matter
        # surface species grouped by site master (Hfo_w, Hfo_s, Hmo_x, Hmo_y, Clay_w, Clay_s, ...)
        site_of = [(re.match(r'^([A-Z][a-z]+_[a-z])', n) or [None, None])[1] for n in s.surf.names]
        site_rows = {g: np.array([x == g for x in site_of], dtype=bool)
                     for g in sorted(set(x for x in site_of if x and not x.startswith(('Hum_', 'Dom_'))))}
        for k in range(nb):
            if not present[k] or fixed[k]:
                continue
            lab = s.comp_elem[s.basis[k]]
            res.partition[lab] = {
                'dissolved': float(s.aq.nu[:, k] @ n_aq),
                'exchange': float(s.ex.nu[:, k] @ n_ex),
                'sorbed': float(s.surf.nu[:, k] @ n_sf),
                'sorbed_hum': float(s.surf.nu[hum_rows, k] @ n_sf[hum_rows]) if hum_rows.any() else 0.0,
                'dom': float(s.surf.nu[dom_rows, k] @ n_sf[dom_rows]) if dom_rows.any() else 0.0,
                'mineral': float(sum(ph_si[p][0][k] * nm[p] for p in active)),
                **{'site:' + g: float(s.surf.nu[rw, k] @ n_sf[rw]) for g, rw in site_rows.items() if rw.any()},
                'free': float(10 ** la[k] * W / 10 ** lg.get(s.basis[k], 0.0))
                if s.basis[k] in lg else float('nan'),
            }
        for alias, sp in s.aliases.items():
            lab = s.comp_elem.get(sp)
            if lab in res.partition and alias not in res.partition:
                res.partition[alias] = res.partition[lab]
        return res
