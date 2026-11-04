"""Offline (one-way) coupling of the TE chemistry to a finished SCEPTER run.

SCEPTER writes, at each of its output times, depth profiles of pH, water
content, flow, major-ion totals (dissolved + exchanged), pCO2 and the
dissolution rate of every mineral. This module reads them and runs the trace
elements through the same soil column:

  release   metal_i += sum_minerals rate(mineral) * content(metal, mineral) * dt
  transport explicit upwind advection of the dissolved metal with SCEPTER's
            water flux (sub-stepped to the water Courant limit)
  chemistry re-partition every cell with solver.Speciator at SCEPTER's pH,
            with the major cations on the TE exchanger and surfaces set in
            equilibrium with SCEPTER's porewater (exchange, Fe-oxide and
            Mn-oxide sorption, minerals)

The host state is held constant over each SCEPTER output interval, so use a
fine output (nstep) for long runs. Units are per m3 of bulk soil: moles, and
water mass W = poro * sat * 1000 kg.

Not yet included: tillage mixing of sorbed metal, dispersion/diffusion,
SOM binding, plant uptake, ageing.
"""
import glob
import os
from dataclasses import dataclass, field

import numpy as np

from . import clay, humic, solver

MW_INRT = 258.162          # SCEPTER 'inrt' is kaolinite (mwtka)
SCEPTER_AQ = {'ca': 'Ca', 'mg': 'Mg', 'k': 'K', 'na': 'Na', 'al': 'Al', 'si': 'Si',
              'no3': 'N(5)'}   # nitrate from N fertiliser (anion, carried if the system has N(5))
MW_EL = {'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38, 'Pb': 207.2, 'Cd': 112.41,
         'Co': 58.933, 'Cr': 51.996, 'As': 74.922, 'Mn': 54.938}


# Site classes for the extraction-definition test. Hmo_x is Tonkin's lower-affinity
# Mn-oxide site (Pb log K far below Hmo_y), grouped with the weak sites.
WEAK_SITES = ['Hfo_w', 'Hmo_x', 'Clay_w']
STRONG_SITES = ['Hfo_s', 'Hmo_y', 'Clay_s']

def _table(path):
    with open(path) as f:
        head = f.readline().split()
        data = np.array([[float(x) for x in l.split()] for l in f if l.strip()])
    return {h: data[:, j] for j, h in enumerate(head)}


@dataclass
class Snapshot:
    time: float
    z: np.ndarray
    dz: np.ndarray
    poro: np.ndarray
    sat: np.ndarray
    q: np.ndarray                 # Darcy flux, m/yr (= pore velocity * poro * sat)
    pH: np.ndarray
    pco2: np.ndarray              # atm
    majors: dict                  # label -> mol/m3 bulk (dissolved + exchanged)
    cec: np.ndarray               # effective CEC (base cations + Al), mol charge / m3 bulk
    rate: dict                    # mineral -> mol/m3/yr (positive = dissolution)
    mw: dict                      # mineral -> g/mol
    cec_total: np.ndarray = None  # SCEPTER's total CEC, mol charge / m3 bulk
    diss: dict = None             # SCEPTER dissolved majors, mol/L

    @property
    def W(self):
        return self.poro * self.sat * 1000.0


def read_run(rundir):
    """All SCEPTER output snapshots of a run, in time order."""
    prof = os.path.join(rundir, 'prof')
    idx = sorted(int(os.path.basename(p)[5:8]) for p in glob.glob(os.path.join(prof, 'rate-[0-9][0-9][0-9].txt')))
    out = []
    for k in idx:
        s = f'{k:03d}.txt'
        bsd = _table(os.path.join(prof, 'bsd-' + s))
        aq = _table(os.path.join(prof, 'prof_aq(tot)-' + s))
        ex = _table(os.path.join(prof, 'prof_ex(tot)-' + s))
        gas = _table(os.path.join(prof, 'prof_gas-' + s))
        rate = _table(os.path.join(prof, 'rate-' + s))
        sld = _table(os.path.join(prof, 'prof_sld-' + s))
        wt = _table(os.path.join(prof, 'prof_sld(wt%)-' + s))
        z = bsd['z']
        edges = np.concatenate([[0.0], 0.5 * (z[1:] + z[:-1]), [z[-1] + 0.5 * (z[-1] - z[-2])]])
        dz = np.diff(edges)
        dens = bsd['dens[g/cm3]']
        poro, sat = bsd['poro'], bsd['sat']
        soil_kg = dens * 1000.0 * (1.0 - poro)
        majors = {lab: ex[sp] for sp, lab in SCEPTER_AQ.items() if sp in ex}
        diss = {lab: aq[sp] for sp, lab in SCEPTER_AQ.items() if sp in aq}
        # Effective CEC for the TE exchanger: the charge SCEPTER holds as
        # exchangeable base cations and Al. Sites SCEPTER fills with H+
        # (its variable-charge, high-H-affinity exchanger) are not offered
        # to the metals, so the two exchange models agree on the majors.
        zch = {'Ca': 2, 'Mg': 2, 'K': 1, 'Na': 1, 'Al': 3}
        cec_eff = sum(zch[l] * np.maximum(majors[l] - poro * sat * 1000.0 * diss[l], 0.0)
                      for l in zch if l in majors)
        # mineral formula weights from mol and wt% profiles, referenced to inrt
        mw = {}
        for m in sld:
            if m in ('z', 'time'):
                continue
            j = int(np.argmax(sld[m]))
            if sld[m][j] <= 1e-12 or sld['inrt'][j] <= 0:
                continue
            ref = wt['inrt'][j] / sld['inrt'][j]
            mw[m] = MW_INRT * (wt[m][j] / sld[m][j]) / ref
        out.append(Snapshot(time=float(bsd['time'][0]), z=z, dz=dz, poro=poro, sat=sat,
                            q=bsd['v[m/yr]'] * poro * sat, pH=aq['aph'],
                            pco2=gas['pco2'], majors=majors,
                            cec=cec_eff, cec_total=bsd['cec[cmol/kg]'] * 0.01 * soil_kg,
                            diss=diss,
                            rate={m: v for m, v in rate.items() if m not in ('z', 'time')},
                            mw=mw))
    return out


@dataclass
class TEResult:
    times: list = field(default_factory=list)
    z: np.ndarray = None
    pools: dict = field(default_factory=dict)   # pool -> list of (cell, metal) arrays, mol/m3
    released: list = field(default_factory=list)  # cumulative mol/m2 per metal
    leached: list = field(default_factory=list)   # cumulative mol/m2 per metal
    metals: list = field(default_factory=list)
    host: list = field(default_factory=list)      # per snapshot: per-cell host totals (mol/m3)
    surfaces: object = None                       # cell -> {name: Surface}
    cec_x: object = None                          # snap -> exchanger charge per cell
    labels: dict = field(default_factory=dict)    # metal -> solver label
    la: list = field(default_factory=list)        # per snapshot: per-cell {basis species: log activity}

    def mgkg(self, pool, k, snap):
        """Pool at output k in mg per kg soil, (cell, metal)."""
        soil_kg = (1 - snap.poro) * 1000 * _dens(snap)
        mw = np.array([MW_EL[m] for m in self.metals])
        return self.pools[pool][k] * mw[None, :] * 1e3 / soil_kg[:, None]


def _dens(snap):
    return getattr(snap, '_dens', 2.65)


def run(sv, snaps, metals, content, *, background=None, fe_mol_per_kg=0.0,
        fe_reactive=0.1, hmo_g_per_kg=0.0, hum_g_per_kg=0.0, extra_totals=None, minerals=(),
        valence_labels=None, courant=0.9, T=None, max_sub=None, cec_mode='total',
        slow=None, implicit_sub=None, dom=None, clay_g_per_kg=0.0):
    """Run trace metals through SCEPTER snapshots.

    sv        : solver.Speciator built with the majors, metals, X, Hfo_w/Hfo_s,
                Hmo_x/Hmo_y as needed, and CO2(g) in phases
    content   : {metal: {mineral: ppm by mass in that mineral}}
    background: {metal: mol/m3 reactive metal already in the soil} (uniform) or None
    fe_mol_per_kg, fe_reactive : soil Fe oxide (mol Fe/kg soil) and its
                ferrihydrite-like fraction, giving D&M Hfo sites
    hmo_g_per_kg : hydrous Mn oxide (Tonkin 2004) per kg soil
    hum_g_per_kg : active humic acid (simplified WHAM VII, chem.humic) per kg
                soil; scalar or one value per cell (e.g. 0.5 x SOM)
    valence_labels: {metal: label passed to the solver}, e.g. {'As': 'As(5)'}
    implicit_sub: if set, transport uses implicit upwind advection of the
                retarded metal (dissolved fraction lagged from the start of each
                sub-step; trace metals are in the linear regime) with this many
                sub-steps per snapshot interval, instead of the explicit Courant-
                limited scheme. For century-scale runs.
    slow      : optional slow (aged / strongly held) pool exchanging with the
                surface-sorbed metal, dS/dt = k_age*sorbed - k_rel*S, e.g.
                {'total': {metal: mol/m3 slow pool at t0}, 'k_age': 1.0}; k_age may be a
                {metal: 1/yr} dict (missing metals do not age). k_rel is set
                per cell and metal at t0 so that the initial state is steady
                (release from S balances ageing plus net leaching from the cell),
                i.e. the control keeps its labile metal under the initial forcing.
    clay_g_per_kg : 2:1 clay (smectite + illite) per kg soil, giving 2SPNE edge
                sites (chem.clay; clay.add_clay must have been applied to the
                database and Clay_s/Clay_w given to the speciator)
    dom       : optional dissolved organic matter, {'doc_ref': mg C/L, 'pH_ref': pH,
                'slope': d log10(DOC)/d pH}. A mobile copy of the humic surface
                (humic.add_dom must have been applied to the database and the Dom_*
                masters given to the speciator) sized from DOC(pH) in each cell's
                porewater. Metal bound to it ('dom' pool) is carried with the water.
    cec_mode  : 'total' offers SCEPTER's whole CEC to the TE exchanger (sites
                SCEPTER leaves H-occupied become HX); 'eff' offers only the
                charge SCEPTER holds as exchangeable base cations and Al.
    """
    def cec_x(snap):
        return snap.cec_total if cec_mode == 'total' else np.maximum(snap.cec, 1e-6)

    lab = {m: (valence_labels or {}).get(m, m) for m in metals}
    s0 = snaps[0]
    nz = len(s0.z)
    soil_kg = (1 - s0.poro) * 1000.0 * 2.65   # only for site densities
    hum = np.broadcast_to(np.asarray(hum_g_per_kg, float), (nz,))
    tot = np.zeros((nz, len(metals)))
    if background:
        for j, m in enumerate(metals):
            tot[:, j] = background.get(m, 0.0)
    res = TEResult(z=s0.z, metals=list(metals),
                   pools={p: [] for p in ('dissolved', 'exchange', 'sorbed', 'sorbed_hum', 'dom', 'mineral', 'total', 'slow',
                                          'sorbed_weak', 'sorbed_strong')})
    res.occ = []          # per snapshot: {site: {element: (nz,) mol on site}, site+'_tot': (nz,) site moles}
    S = np.zeros((nz, len(metals))); k_rel = None
    if slow:
        for j, m in enumerate(metals):
            S[:, j] = slow['total'].get(m, 0.0)
        k_age = slow.get('k_age', 1.0)
        if isinstance(k_age, dict):       # per-metal ageing rate; metals left out do not age
            k_age = np.array([float(k_age.get(m, 0.0)) for m in metals])
    cum_rel = np.zeros(len(metals))
    cum_lea = np.zeros(len(metals))
    guesses = [None] * nz
    t_prev = 0.0
    diss = None

    def surfaces(i, snap=None):
        surf = {}
        if dom and snap is not None:
            doc = humic.doc_mg_per_L(float(snap.pH[i]), dom['doc_ref'], dom['pH_ref'], dom['slope'])
            surf['Dom'] = humic.dom_surface(humic.OM_PER_DOC * doc * 1e-3 * float(snap.W[i]))
        fe = fe_mol_per_kg * fe_reactive * soil_kg[i]
        if fe > 0:
            surf['Hfo'] = solver.Surface({'Hfo_w': 0.2 * fe, 'Hfo_s': 0.005 * fe}, 600.0, 89.0 * fe)
        if hmo_g_per_kg > 0:
            g = hmo_g_per_kg * soil_kg[i]
            surf['Hmo'] = solver.Surface({'Hmo_x': g * 2.1e-3 * 0.64, 'Hmo_y': g * 2.1e-3 * 0.36}, 746.0, g)
        if clay_g_per_kg > 0:
            surf['Clay'] = clay.surface(clay_g_per_kg * soil_kg[i])
        h = hum[i] * soil_kg[i]
        if h > 0:
            surf['Hum'] = humic.surface(h)
        return surf

    def solve(snap, i, totals, W):
        return sv.solve(totals, pH=float(snap.pH[i]), W=W, T=T or 288.15,
                        cec={'X': float(cec_x(snap)[i])}, surfaces=surfaces(i),
                        fix_si=[('CO2(g)', float(np.log10(max(snap.pco2[i], 1e-6))), 'CO3-2')],
                        minerals=list(minerals))

    known = set(sv.s.comp_elem.values()) | set(sv.s.aliases)

    def host_totals(snap):
        """Major-element totals per cell such that, with the TE exchanger and
        surfaces, the porewater has SCEPTER's dissolved majors. Solved once per
        snapshot by equilibrating the exchanger and surfaces with SCEPTER's
        porewater in a huge water mass (as PHREEQC -equilibrate does)."""
        big = 1e6
        out = []
        for i in range(nz):
            sol = {k: max(v[i], 1e-20) for k, v in snap.diss.items() if k in known}
            if extra_totals:
                sol.update(extra_totals)
            r = solve(snap, i, {k: v * big for k, v in sol.items()}, big)
            t = {}
            for k, v in sol.items():
                p = r.partition[k]
                t[k] = v * snap.W[i] + p['exchange'] + p['sorbed']
            out.append(t)
        return out

    def equilibrate(snap, hostt):
        part = {p: np.zeros_like(tot) for p in ('dissolved', 'exchange', 'sorbed', 'sorbed_hum', 'dom', 'mineral',
                                                 'sorbed_weak', 'sorbed_strong')}
        occ = {}
        for i in range(nz):
            totals = dict(hostt[i])
            for j, m in enumerate(metals):
                totals[lab[m]] = max(tot[i, j], 1e-30)
            r = sv.solve(totals, pH=float(snap.pH[i]), W=float(snap.W[i]), T=T or 288.15,
                         cec={'X': float(cec_x(snap)[i])}, surfaces=surfaces(i, snap),
                         fix_si=[('CO2(g)', float(np.log10(max(snap.pco2[i], 1e-6))), 'CO3-2')],
                         minerals=list(minerals), guess=guesses[i])
            guesses[i] = dict(r.la)
            for j, m in enumerate(metals):
                p = r.partition[lab[m]]
                for k in part:
                    part[k][i, j] = p.get(k, 0.0)
                part['sorbed_weak'][i, j] = sum(p.get('site:' + g, 0.0) for g in WEAK_SITES)
                part['sorbed_strong'][i, j] = sum(p.get('site:' + g, 0.0) for g in STRONG_SITES)
            # site occupancy by every element (metals, Ca, Mg, Al, ...), and site totals
            for sf in surfaces(i, snap).values():
                for g, ntot in sf.sites.items():
                    if g not in WEAK_SITES + STRONG_SITES:
                        continue
                    d = occ.setdefault(g, {})
                    d.setdefault('_tot', np.zeros(nz))[i] = ntot
                    for el, p in r.partition.items():
                        v = p.get('site:' + g, 0.0)
                        if v > 0:
                            d.setdefault(el, np.zeros(nz))[i] = v
        equilibrate.occ = occ
        part['sorbed'] -= part['dom']          # 'sorbed' = solid surfaces only
        return part

    for snap in snaps:
        dt_tot = snap.time - t_prev
        if dt_tot <= 0:
            continue
        # water Courant limit on the explicit advection
        dt_c = courant * np.min(snap.W / 1000.0 * snap.dz / np.maximum(snap.q, 1e-12))
        nsub = max(1, int(np.ceil(dt_tot / dt_c))) if not implicit_sub else implicit_sub
        if max_sub:
            nsub = min(nsub, max_sub)
        dt = dt_tot / nsub
        src = np.zeros_like(tot)
        for j, m in enumerate(metals):
            for mineral, ppm in content.get(m, {}).items():
                r = snap.rate.get(mineral)
                if r is None or mineral not in snap.mw:
                    continue
                # mol metal per mol mineral
                f = ppm * 1e-6 * snap.mw[mineral] / MW_EL[m]
                src[:, j] += np.maximum(r, 0.0) * f
        hostt = host_totals(snap)
        for _ in range(nsub):
            tot += src * dt
            cum_rel += (src * snap.dz[:, None]).sum(0) * dt
            part = equilibrate(snap, hostt)
            mob = part['dissolved'] + part['dom']             # DOM-bound metal moves with the water
            c = mob / snap.W[:, None]                         # mol/kgw
            flux = snap.q[:, None] * 1000.0 * c               # mol/m2/yr leaving each cell downward
            inflow = np.vstack([np.zeros((1, len(metals))), flux[:-1]])
            if slow:
                if k_rel is None:     # steady initial state
                    net_out = np.maximum((flux - inflow) / snap.dz[:, None], 0.0)
                    k_rel = (k_age * part['sorbed'] + net_out) / np.maximum(S, 1e-30)
                    k_rel = np.where(S > 0, k_rel, 0.0)
                    res.k_rel = k_rel.copy()
                age = part['sorbed'] * (1 - np.exp(-k_age * dt))
                rel = S * (1 - np.exp(-k_rel * dt))
                S += age - rel
                tot += rel - age
            if implicit_sub:
                # backward-Euler upwind on the total, c = fd * T / W (fd lagged)
                fd = mob / np.maximum(tot, 1e-300)
                a = snap.q[:, None] * 1000.0 * fd / snap.W[:, None] * dt / snap.dz[:, None]
                newt = np.zeros_like(tot)
                for i in range(nz):
                    up = a[i - 1] * newt[i - 1] * snap.dz[i - 1] / snap.dz[i] if i > 0 else 0.0
                    newt[i] = (tot[i] + up) / (1.0 + a[i])
                out_flux = a[-1] * newt[-1] * snap.dz[-1] / dt
                cum_lea += out_flux * dt
                tot[:] = newt
                continue
            tot += (inflow - flux) * dt / snap.dz[:, None]
            cum_lea += flux[-1] * dt
        part = equilibrate(snap, hostt)
        t_prev = snap.time
        res.times.append(snap.time)
        for k in part:
            res.pools[k].append(part[k])
        res.pools['total'].append(tot.copy())
        res.occ.append(getattr(equilibrate, 'occ', {}))
        res.pools['slow'].append(S.copy())
        res.released.append(cum_rel.copy())
        res.leached.append(cum_lea.copy())
        res.host.append(hostt)
        res.la.append([dict(g) if g else {} for g in guesses])   # log activities of basis species per cell
    res.surfaces = surfaces
    res.cec_x = cec_x
    res.labels = lab
    return res
