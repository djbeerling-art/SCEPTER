"""One-way coupled trace-element column (SCEPTER-TE Phase 1 prototype).

The host model (SCEPTER, or PHREEQC in the benchmark) supplies, for every
cell and step, the converged pH, water mass and system totals of the major
elements (dissolved + exchangeable). The TE module carries only the trace
metals: each step it advects the dissolved part of every metal one cell
down (pure advection, Courant number 1, as PHREEQC TRANSPORT shifts), adds
the inflow, and re-equilibrates every cell with the host state. Metals are
too dilute to change pH or the majors, which is what makes one-way coupling
valid; the benchmark tests that assumption as well as the chemistry.
"""
from dataclasses import dataclass, field

import numpy as np



@dataclass
class ColumnResult:
    metals: list
    pools: dict = field(default_factory=dict)   # pool -> array (step, cell, metal)
    effluent: list = field(default_factory=list)  # (step, {metal: moles leaving})


def run(sv, host, metals, init_totals, inflow_moles, *, cec=None, surfaces=None,
        fix_si=(), minerals=(), T=298.15, pe=None):
    """sv: solver.Speciator; host[k][i] = dict(pH, W, totals={label: mol},
    optional pe) for
    step k (k = 0 is the initial equilibrium) and cell i; init_totals[i] =
    {metal: mol}; inflow_moles = {metal: mol per shift entering cell 0}."""
    nstep, ncell = len(host), len(host[0])
    out = ColumnResult(metals=list(metals))
    pools = {p: np.zeros((nstep, ncell, len(metals)))
             for p in ('dissolved', 'exchange', 'sorbed', 'mineral', 'total')}
    tot = [dict(t) for t in init_totals]
    diss = [None] * ncell
    guesses = [None] * ncell
    for k in range(nstep):
        if k > 0:
            moving = [{m: diss[i][m] for m in metals} for i in range(ncell)]
            out.effluent.append((k, moving[-1]))
            for i in range(ncell):
                src = inflow_moles if i == 0 else moving[i - 1]
                for m in metals:
                    tot[i][m] = tot[i][m] - moving[i][m] + src.get(m, 0.0)
        for i in range(ncell):
            h = host[k][i]
            totals = dict(h['totals'])
            totals.update({m: max(tot[i][m], 1e-30) for m in metals})
            r = sv.solve(totals, pH=h['pH'], W=h['W'], T=T, pe=h.get('pe', pe), cec=cec,
                         surfaces=surfaces, fix_si=fix_si, minerals=minerals,
                         guess=guesses[i])
            guesses[i] = {b: v for b, v in r.la.items()}
            diss[i] = {}
            for j, m in enumerate(metals):
                p = r.partition[m]
                diss[i][m] = p['dissolved']
                for pool in ('dissolved', 'exchange', 'sorbed', 'mineral'):
                    pools[pool][k, i, j] = p[pool]
                pools['total'][k, i, j] = tot[i][m]
    out.pools = pools
    return out
