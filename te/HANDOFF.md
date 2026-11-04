# SCEPTER-TE handoff (updated 3 Nov 2026)

Goal (David): short paper on trace-element (Ni, Cr, Cu, Zn, Cd, Pb) build-up under basalt
enhanced weathering in AERATED upland soils, as a mechanistic answer to Dupla et al. 2023
(EJSS 74:e13343). 3 temperate + 3 tropical soils, several basalts, 10-40 t/ha/yr, 10-40 yr.

Plan doc: https://claude.ai/artifact/4V212c94Do8WZqrCrWxDZ9

Design: one-way coupled TE module inside upstream (cropland) SCEPTER, called after each
converged time step (metals too dilute to feed back on pH/majors). Processes: host-mineral
release, aqueous speciation, exchange, DLM surface complexation on Fe oxides, humic binding,
solubility caps, Cr redox (Mn oxides), plant uptake. Prototype chemistry in Python vs PHREEQC first.

Status:
- Build: apt-get update && apt-get install gfortran libopenblas-dev; make --file=makefile_gmd24 -> scepter_gmd24.
  pip install scipy pandas matplotlib phreeqpython.
- te/upland.py: spin-up (inert + g2 OM + CEC, Ca rain sets pH; 10 kyr, w=1e-3, ~25 s) and
  amend() (restart, Hillhouse la/agt/fo/anl as dust, tilling mix, 2 yr ~40 s).
  Porewater pH vs Ca rain (10 C, q 0.25 m/yr, CEC 12, no N): 5e-6 -> 5.5, 1e-4 -> 6.4, 5e-4 -> 7.4.
  N fertiliser (NH4NO3) without base drives pH to 2.4: needs liming or skip N.
- Gotchas: rstrt must be an absolute path; every aq species needs a nonzero rain conc
  (1e-20; Fe3+ 1e-21) or the run goes NaN; cec.in needs 5 logK values + beta.

Scope decisions (David, 3 Nov; also in the plan doc):
- Staged build. Phase 1 is Ni + Cu (release, speciation, Fe-oxide sorption, exchange, transport),
  validated against PHREEQC, Kelland 2020 and Vienne 2022. Phase 2 is Zn, Cd, Pb + OM binding,
  DOC transport and plant uptake. Phase 3 is Cr (Cr(III) release, Mn-oxide and redox-kinetics
  sensitivity of Cr(VI)).
- Labile = exchangeable + weakly/reversibly sorbed. Strong-site, high-affinity SOM and aged metal
  are reported separately; the boundary between them is a sensitivity parameter.
- Headline result: whether Dupla-style screening and the mechanistic model agree on the risk class
  for each soil x basalt x dose (not threshold dates).
- Dupla PDF, plan .docx and Beerling et al. 2024 PNAS (Corn Belt field trial, T2 site, 50 t/ha/yr x 4 yr,
  trace metals in SI Figs S16-S18 + Dataset S1, not yet obtained) are in /mnt/project-files/te/refs/.

Chemistry prototype (te/chem, see te/chem/README.md):
- Reads minteq.v4.dat + te/chem/data/te_extra.dat directly, so its constants are PHREEQC's.
  Fixed-pH Newton solver: aqueous speciation, Gaines-Thomas exchange, D&M DLM on Hfo (w/s),
  fixed-SI components (CO2(g)) and mineral solubility caps (active set). Redox is 'decoupled'
  by default (Cr(3)/Cr(6) as separate components) or 'coupled' at a fixed pe. ~10 ms/solve.
- python3 te/tests/validate_phreeqc.py runs batch pH edges (pH 4-8) against PHREEQC:
  Ni/Cu/Zn/Cd/Pb/Cr(VI)/Cr(III) on Hfo, 0.5 kg water, a soil-like mix with exchange and
  pCO2 at 10 and 25 C, and malachite/cerussite/smithsonite caps. Worst pool-fraction error
  1.2e-5, free ion within 2e-4. Results: te/tests/out/validation.md, pH_edges.png.
- PHREEQC gotchas: MOL() returns mol/kgw for exchange and surface species too. minteq's
  "-gamma 0 0" means the DH limiting law (log gamma 0 for neutral species), not Davies.
  Don't fix a high pe with Fix_pe/O2(g) (PHREEQC fails); use O2(g) SI instead.
  Cu needs Cu(2) as its decoupled valence.

Field Ni check (te/studies/ni_batch_field.py, output alongside), closed batch, equilibrium only:
- Corn Belt (T2): pH +0.2-0.4 with higher Ca gives NH4OAc-like Ni 0.5-0.7x control at 10-30%
  reactive Fe; observed ~0.55x (Fig S18). Equilibrium sorption onto Fe oxide explains it.
- Sabah (P2): observed 0.23x (2018->2024, pH 5.5->5.74). Closed batch gives only 0.7-0.95x, but
  higher Ca triples porewater Ni (~100 -> 300 ug/L). At ~1 m/yr drainage that leaches ~3 kg Ni/ha/yr,
  about the whole 0-10 cm exchangeable pool, so Ca displacement + leaching is the likely mechanism.
  Tested in the transport column (below). Al exchange and SOM binding are not in yet.

Transport column (done 3 Nov): `te/chem/column.py` is the one-way coupled column. The host gives
pH, water and major-element system totals per cell and step; the TE module advects the dissolved
metal (Courant 1, as PHREEQC shifts) and re-partitions it. `te/tests/validate_column.py` runs
an acid soil column (20 cells, exchanger + Hfo, Ni and Cu loaded) flushed with Ca-Mg-HCO3 water
in PHREEQC TRANSPORT and in the Python column driven by PHREEQC's host state: worst error 5e-5,
so one-way coupling holds for trace Ni and Cu. Note that the Ca front pushes Ni ahead of it
(effluent Ni peaks at ~3.7x the initial porewater).
Sabah pH-versus-leaching runs live in the project files (unpublished data, not in the repo).
They point to Al exchange, Mn oxides and SOM as missing pieces in acid soils.

Al exchange (done 3 Nov): AlX3 and AlOHX2 (phreeqc.dat constants) are in te_extra.dat, and an
acid_soil_Al case with gibbsite checks them. Gotcha: PHREEQC gives an exchange species with
`-gamma` the activity coefficient of the LAST aqueous token of its reaction rewritten in master
species, scaled by that token's coefficient, and uses Davies when a = b = 0. So AlOHX2 gets the
H+ term with coef -1 and a positive log g (system._ex_gamma). The column test now carries Ni, Cu,
Pb and Cd with Al exchange and gibbsite; worst error 5e-5.
Field check (acid soil, project files only): gibbsite-controlled Al matches the low NH4OAc Al;
Fe-oxide sorption + exchange reproduce the Pb fall and flat Cd, but not the Ni fall (Ni tracks Mn
and Co) or flat Cu, which points to Mn oxides (Ni) and SOM (Cu).

Next steps:
1. Arsenic batch chemistry done (edge_AsV, edge_AsIII match PHREEQC). Next for As: competition
   from phosphate, silicate and bicarbonate on Hfo (basalt releases Si), then As in the column.
   Mn oxide done: Tonkin et al. 2004 HMO DLM (Hmo_x/Hmo_y) for Ni, Co, Cu, Zn, Cd, Pb, Mn,
   validated against PHREEQC. Gotchas: minteq's Co master is Co+3 and Mn has several valences,
   so carry Co(2)/Mn(2) as decoupled valences. Equilibrium HMO reproduces the Sabah-like Co and
   Mn falls but only part of the Ni fall, and every oxide model over-predicts the Cu fall:
   next pieces are SOM binding (Cu, Pb) and a slow Ni ageing/incorporation pool on Mn oxide.
2. Offline coupling done (te/chem/scepter_te.py): reads a SCEPTER run, releases metals from
   dissolving host minerals, advects and re-partitions them; mass balance closes. First Vermont
   run (study in project files) gives labile Ni 0.9x, Zn 0.8x, Cd 1.0x after 3 yr, with Cu and Pb
   falling too far (no SOM binding yet), and pH rising ~0.4 vs the observed 0.15-0.24.
   SCEPTER control runs with dose 0 stall in the PSD solver: use the spin-up state instead.
   With site-matched feedstocks (upland.py ROCKS/PSD 'pioneer', 'blueridge'; studies in project
   files): Vermont 3 yr gives pH +0.1 and all labile metals 0.9-1.0x (field: no change). Corn
   Belt 4 x 50 t/ha gives labile Ni 0.51x (field ~0.55x), Cd and Pb falling (field: fall), but Cu
   0.15x (field: no change -> SOM) and pH +1.1 (no N fertiliser acidity in the spin-up).
   Not yet in: the feedstock's own exchangeable Ca (Pioneer Valley ECEC 50 meq/100 g).
3. Per-soil Hfo site densities from Fe-oxide content: ~0.6 mol Fe/kgw if 1.5% Fe oxide were all ferrihydrite at
   field moisture, far above the 1 mmol used in the tests, so sorption will dominate over
   exchange in most soils.
4. Fortran port into SCEPTER as the one-way module.
