# te/chem: trace-element equilibrium chemistry (Python prototype)

This is the chemistry core of the SCEPTER-TE module, prototyped in Python and
benchmarked against PHREEQC before it is ported to Fortran. At each SCEPTER
time step and depth it takes the converged local state (pH, T, water mass,
major-ion totals, pCO2, CEC, Fe-oxide sites) and the trace-metal totals. It
returns how each metal splits between dissolved (free ion and complexes),
exchangeable, sorbed on Fe oxide and precipitated.

| file | role |
|---|---|
| `thermo.py` | Reads PHREEQC-format databases (master species, aqueous, exchange and surface species, phases, log K, delta H, analytic and -gamma). |
| `system.py` | Reduces every species to a basis of components. In `redox='decoupled'` mode (the SCEPTER-TE default) chosen valence states, such as Cr(3) and Cr(6), are separate components, so redox can be kinetic. `'coupled'` links them through a fixed pe, as PHREEQC does. |
| `solver.py` | Newton solver at fixed pH (and fixed pe in coupled mode), with an analytic Jacobian. It covers aqueous speciation (Davies or WATEQ Debye-Hückel, as PHREEQC), Gaines-Thomas exchange, a Dzombak & Morel diffuse-layer model on Hfo strong and weak sites, components fixed by saturation index (such as CO2(g)), and solubility caps through an active set of minerals. About 10 ms per solve for the full 12-element soil system. |
| `column.py` | One-way coupled 1-D column: the host (SCEPTER, or PHREEQC in the test) gives pH, water and major totals per cell; the module advects dissolved metal and re-partitions it each step. |
| `phreeqc_ref.py` | Runs the same batch problem (`run`) or advective column (`run_transport`) in PHREEQC (IPhreeqc via phreeqpython) and returns PHREEQC's system totals, so both codes get identical inputs. |
| `data/minteq.v4.dat` | USGS PHREEQC database, used unchanged. It provides the D&M Hfo constants for Ni, Cu, Zn, Cd, Pb, Cr(III) and chromate. |
| `data/te_extra.dat` | Exchange constants (Appelo & Postma 2005) and the Fix_H+/Fix_pe helper phases. Both PHREEQC and the prototype read this file. |

Validation: `python3 te/tests/validate_phreeqc.py` (needs `pip install phreeqpython`).
It writes `te/tests/out/validation.md` and `pH_edges.png`, and exits non-zero
if any pool fraction differs from PHREEQC by more than 0.01.
`python3 te/tests/validate_column.py` checks the column against PHREEQC
TRANSPORT and writes `te/tests/out/column_validation.md` and `column_profiles.png`.

Conventions match PHREEQC. Aqueous activities are molal and surface master
"activities" are moles. Neutral species with `-gamma` get log γ = bμ; neutral
species without it get 0.1μ. Charged species without `-gamma` use Davies.
Exchange species get the WATEQ γ only when `-gamma` is given. The dielectric
constant is Bradley & Pitzer (1979).

Not yet included: humic/SOM binding and DOC (Phase 2), goethite sites, Cr
redox kinetics (Phase 3), and H+/Al3+ exchange in strongly acid soils (HX is
in te_extra.dat, Al is not).
