# Riceland SCEPTER: process specification (v15)

This document lists every process added to SCEPTER for flooded rice paddies: rate laws, stoichiometry,
parameters, units, temperature and pH treatment, boundary schedules and the calibration behind each
value. It is written so the processes can be re-implemented in another reactive-transport code (for
example as PHREEQC `RATES`/`KINETICS` blocks) without porting the Fortran. Code references are to
`scepter.f90` and `paddy/run_flooded_column.py` on branch `paddy/fe3-reduction` (model version v15,
calibration file `paddy/krxn_calibrated_v15.in`).

Notation: "Monod(X, K)" = X/(X + K); "Inh(X, Ki)" = Ki/(X + Ki); φ porosity; S_w water saturation;
K_H Henry's constant (mol L-1 atm-1); [X] aqueous concentration (M).

---

## 1. Units and conversion to PHREEQC

### 1.1 SCEPTER units

| Quantity | SCEPTER unit | Notes |
|---|---|---|
| Time | yr (365 d = 3.1536e7 s, `sec2yr`) | |
| Depth | m | |
| Solids (OM `g2`, minerals, urea) | mol m-3 of bulk soil | |
| Aqueous species | mol L-1 of pore water (M); storage term φ·S_w·1000·[X] mol m-3 bulk | `scepter.f90` aqueous residual |
| Gases | partial pressure, atm; total = air-filled + dissolved (K_H·p) | at S_w = 1 all of it is dissolved |
| Reaction rates `rxnext` | mol of reaction extent per m3 bulk per yr | enter every mass balance as stoich × rate |
| Root/plant fluxes | mol m-2 yr-1 of field, spread over the root profile | see §5 |
| OM input | g C m-2 yr-1 (default 450) into the top 0.15 m | |

For the OM-consuming reactions the extent is **1 mol of CH2O consumed**; every stoichiometric
coefficient below is per mol CH2O, so `vmax` is a rate of OM loss.

### 1.2 Conversion factors

For a cell of bulk volume 1 m3:

* pore water: W = φ·S_w·1000 L ≈ φ·S_w·1000 kg water (density taken as 1 kg/L)
* rate: **R_PHREEQC (mol kgw-1 s-1) = R_SCEPTER (mol m-3 yr-1) / (φ·S_w·1000 · 3.1536e7)**
* solid or kinetic-reactant amount: **M_PHREEQC (mol kgw-1) = M_SCEPTER (mol m-3) / (φ·S_w·1000)**
* first-order constants (1/yr) divide by 3.1536e7 to give 1/s
* aqueous M ≈ mol/kgw for these dilute pore waters

At the reference setting (φ = 0.5, flooded S_w = 1) the factor is 1/(500 · 3.1536e7) = 6.34e-11.
During drains and the fallow S_w < 1 in the upper cells, so the factor changes with depth and time;
compute it per cell. Rate laws written in terms of solid amounts (e.g. Monod on g2 with K = 100 mol m-3)
must have the half-saturation constant converted with the same factor, or be evaluated in SCEPTER units
inside the RATES block and the result converted.

---

## 2. Species

Types: aq = aqueous (SCEPTER tracks totals; speciation is by equilibrium), s = solid,
g = gas (SCEPTER "gas" species are total dissolved + gaseous). "phreeqc.dat" refers to the database
shipped with PHREEQC 3; llnl.dat may differ.

| SCEPTER name | Meaning | Type | In phreeqc.dat? | Notes for PHREEQC |
|---|---|---|---|---|
| `g2` | reactive OM, CH2O with N:C = n2c (0.03) | s | no | kinetic reactant, formula CH2O + NH3 at N:C 0.03 |
| `urea` | urea fertiliser | s | no | kinetic reactant CO(NH2)2, hydrolyses in days |
| `amfe3` | ferrihydrite | s | yes, `Fe(OH)3(a)` | |
| `gt` | goethite | s | yes, `Goethite` | no precipitation allowed in the paddy runs |
| `gr` | green rust (Fe3(OH)8-type, Ponnamperuma 1972) | s | **no** | needs a PHASES entry consistent with §4 (log K −9.98 for Fe2+ + 2 Fe(OH)3 + 2 H2O = Fe3(OH)8 + 2 H+), or treat `fe2gr` as a kinetic reaction using SCEPTER's Ω directly |
| `sid` | siderite | s | yes, `Siderite` | not nucleated in production (seed 0) |
| `fes` | amorphous FeS | s | yes, `FeS(ppt)` / `Mackinawite` | SCEPTER log K 3.915 for Fe2+ + HS- = FeS + H+ |
| `cc` | calcite | s | yes, `Calcite` | nucleation seed 0.001 |
| `inrt` | inert soil matrix carrying the CEC | s | n/a | use EXCHANGE |
| basalt minerals (fo, en/di, ab, an, ka, gb, amsi …) | | s | mostly yes | SCEPTER standard kinetics; rate scaling in §8 |
| `fe2` | total Fe(II) | aq | yes, Fe(+2) | |
| `fe3` | total Fe(III) | aq | yes, Fe(+3) | trace only |
| `ca`, `mg`, `na`, `k`, `cl`, `si`, `al` | | aq | yes | |
| `so4` | total sulfate | aq | yes, S(+6) | |
| `hs` | total sulfide | aq | yes, S(−2) | |
| `no3` | nitrate | aq | yes, N(+5) | |
| `pco2` | total inorganic C (DIC) | g | yes, C(+4), `CO2(g)` | |
| `po2` | O2 | g | yes, `O2(g)`; decoupled `Oxg` also present | |
| `pch4` | CH4 | g | yes, C(−4) and `CH4(g)`; decoupled `Mtg`/`Mtg(g)` also present | use `Mtg` to keep CH4 out of redox equilibrium |
| `pnh3` | total ammonia (NH3 + NH4+) | g | yes, N(−3), `NH3(g)` | |
| `pn2o` | N2O | g | **no** | define a new non-redox master species (e.g. `N2o`) with SOLUTION_SPECIES and a `N2O(g)` PHASE (K_H ≈ 0.024 M/atm at 25 °C) |
| N2 | product of denitrification | — | (`Ntg`) | not tracked in SCEPTER (lost from the column) |

**Redox: important for PHREEQC.** SCEPTER has no redox equilibrium. Every electron transfer
(Fe(III)/Fe(II), SO4/HS, CO2/CH4, NO3/N2O/N2, O2 consumption) is a kinetic reaction from §3–§4. In
PHREEQC, coupled redox states would re-equilibrate instantly through pe and override these rates. Use
the decoupled master species (`Mtg`, `Oxg`, `Ntg`, a new `N2o`) and check that Fe(2)/Fe(3) and
S(6)/S(−2) are not driven to pe equilibrium in the configuration ARTEMIS uses. Check the database
version you use, because the decoupled species were added in later phreeqc.dat releases.

---

## 3. Common rate modifiers

| Modifier | Form | Value | Source |
|---|---|---|---|
| OM limitation | Monod(g2, 100 mol m-3) on all OM-consuming reactions | | SCEPTER default |
| Saturation factor | f_sat = min(2·S_w, 1) on the anaerobic reactions | | this work |
| O2 inhibition (anaerobic pathways) | Inh([O2], 2.5e-7 M) (`--o2inh`) | | v7 calibration (§9) |
| NO3 inhibition (Fe, SO4, CH4 pathways) | Inh([NO3-], 11.3e-5 M) | | Maggi et al. 2008 |
| Plant factor | 1 in season, 0 in fallow (`plant_temp.in`) | | §6 |

---

## 4. Reactions

Stoichiometry is given twice: as SCEPTER writes it (components only, per mol of extent; no H+ or H2O,
see §7) and as a charge-balanced form for PHREEQC.

### 4.1 Carbon, iron and sulfur

| Reaction | SCEPTER stoichiometry (per extent) | Balanced form | Rate law (mol m-3 yr-1) | Parameters (v15) |
|---|---|---|---|---|
| `g2o2` aerobic respiration | −1 g2, −1 O2, +1 DIC, +n2c NH3 | CH2O + O2 → CO2 + H2O | v·Monod(g2)·Monod([O2], 5e-6 M) | v 300 |
| `g2fe3` ferrihydrite reduction | −1 g2, −4 amfe3, +4 Fe(II), +1 DIC, +n2c NH3 | CH2O + 4 Fe(OH)3 + 8 H+ → CO2 + 4 Fe2+ + 11 H2O | v·Monod(g2)·Fe_pool/(Fe_pool + 100)·f_O2·f_NO3·f_sat; Fe_pool = amfe3 + 2·f_gr·gr (mol m-3) | v 300 |
| `g2gr` green-rust reduction | −1 g2, −2 gr, +6 Fe(II), +1 DIC | CH2O + 2 Fe3(OH)8 + 12 H+ → CO2 + 6 Fe2+ + 15 H2O | shares the `g2fe3` expression; GR enters Fe_pool with bioavailability f_gr | f_gr 0.3 |
| `g2gt` goethite reduction | −1 g2, −4 gt, +4 Fe(II), +1 DIC | CH2O + 4 FeOOH + 8 H+ → CO2 + 4 Fe2+ + 7 H2O | v·Monod(g2)·Monod(gt, 100)·f_O2·f_NO3·f_sat | v 10 (crystalline Fe ~30× slower; Roden & Zachara 1996) |
| `g2so4` sulfate reduction | −1 g2, −0.5 SO4, +0.5 HS, +1 DIC | CH2O + ½ SO4 2- + ½ H+ → CO2 + ½ HS- + H2O | v·Monod(g2)·Monod([SO4], 1e-4 M)·f_O2·f_NO3·f_sat | v 300 |
| `g2ch4` methanogenesis | −1 g2, +0.5 CH4, +0.5 DIC | CH2O → ½ CH4 + ½ CO2 | v·Monod(g2)·30/(Fe_pool + 30)·Inh([SO4], 1e-4 M)·f_O2·f_NO3·f_sat | v 1630, Fe(III) Ki 30 mol m-3 |
| `ch4o2` methanotrophy | −1 CH4, −2 O2, +1 DIC | CH4 + 2 O2 → CO2 + 2 H2O | v·Monod([CH4], 5e-6 M)·Monod([O2], 2e-5 M) | v 500 mol m-3 yr-1 |
| `fe2o2` Fe(II) oxidation | −1 Fe(II), −¼ O2, +1 Fe(III)aq | Fe2+ + ¼ O2 + H+ → Fe3+ + ½ H2O (then Fe(OH)3 precipitates kinetically) | (k_SS·([OH-])² + k0)·[Fe2+]·pO2 on dissolved and exchangeable Fe(II) | k_SS 8e13 M-2 min-1 atm-1 (Singer & Stumm 1970), k0 1e-7 min-1 atm-1 |
| `fe2gr` green-rust pH-stat | −1 Fe(II), −2 amfe3, +1 gr | Fe2+ + 2 Fe(OH)3 + 2 H2O → Fe3(OH)8 + 2 H+ | v·f_sat·ln Ω with Ω = 10^(−9.98)·[Fe2+]free/[H+]² (Ω > 1 forms GR); × Monod(amfe3, 1 mol m-3) when forming, × Monod(gr, 1 mol m-3) when dissolving | v 1e4 |
| `groxy` green-rust oxidation | −1 gr, −¼ O2, +3 amfe3 | Fe3(OH)8 + ¼ O2 + ½ H2O → 3 Fe(OH)3 | 1e7·gr·K_H,O2·pO2 | |
| `fe2fs` FeS precipitation | −1 Fe(II), −1 HS, +1 fes | Fe2+ + HS- → FeS + H+ | 1e3·ln Ω, log K 3.915 | |
| `fesox` FeS oxidation | −1 fes, −2 O2, +1 Fe(II), +1 SO4 | FeS + 2 O2 → Fe2+ + SO4 2- | 1e6·fes·[O2] | |
| `sidox` siderite oxidation | −1 sid, −¼ O2, +1 amfe3, +1 DIC | FeCO3 + ¼ O2 + 3/2 H2O → Fe(OH)3 + CO2 | 1e5·sid·[O2] | |

Green rust stands in for the pH buffering of flooded soils near pH 6.5–7 (Ponnamperuma 1972); log K
−9.98 from Bourrié et al. 1999. Goethite and hematite are not allowed to precipitate (ferrihydrite would
otherwise ripen during spin-up).

### 4.2 Nitrogen

| Reaction | SCEPTER stoichiometry | Balanced form | Rate law | Parameters (v15) |
|---|---|---|---|---|
| `amo2o` nitrification | −1 NH3, −2 O2, +1 NO3 | NH4+ + 2 O2 → NO3- + 2 H+ + H2O | v·Monod([NH4], 14e-5 M)·Monod([O2], 2.41e-5 M)·exp(−0.5(pH − 7)²) | v 30 yr-1 (krxn file) |
| `g2n21` denitrification to N2O | −1 g2, −1 NO3, +0.5 N2O, +1 DIC | CH2O + NO3- + H+ → CO2 + ½ N2O + 3/2 H2O | v·Monod(g2)·Monod([NO3], 11.3e-5 M)·f_O2·exp(−0.5(pH − 7)²) × r | v 300 |
| `g2n0` denitrification to N2 | −1 g2, −0.8 NO3, +1 DIC (N2 not tracked) | CH2O + 4/5 NO3- + 4/5 H+ → CO2 + 2/5 N2 + 7/5 H2O | 1.25·(1 − r) × the g2n21 expression without r | |
| N2O yield r | r = 0.1 + 0.9/(1 + 10^(pH − 6.2)) (`--n2oy=6.2:1:0.1`) | | | Simek & Cooper 2002; Liu et al. 2010; Russenes et al. 2016 |
| `g2n2p` N2O reduction | −1 g2, −2 N2O, +1 DIC | CH2O + 2 N2O → CO2 + 2 N2 + H2O | v·Monod(g2)·Monod([N2O], 5e-6 M)·Inh([NO3], 1e-3 M)·f_O2·1/(1 + ([H+]/10^−5.5)²) | v 300, pH50 5.5 |
| `nh4up` plant N uptake | −1 NH3 | NH4+ removed | 3.5 mol N m-2 yr-1 · Monod([NH4], 5e-5 M), root profile | |

Fertiliser: urea, 120 kg N ha-1 per season, added evenly over the 113-day season with the OM input.
Residue N enters with OM at N:C 0.03 (`--n2c`). Nitrification and denitrification forms follow Maggi
et al. 2008.

---

## 5. Plant and root terms

Each flux F (mol m-2 yr-1) is spread over the root profile as F·exp(−z/0.1 m)/∫exp(−z/0.1 m)dz and
multiplied by the plant factor (0 in fallow).

| Term | Effect | Value |
|---|---|---|
| `rtco2` root respiration | + DIC | 60 mol C m-2 yr-1 |
| `rtoxy` radial O2 loss | + O2, × Si feedback f_Si | 3.65 mol O2 m-2 yr-1 (10 mmol m-2 d-1) |
| Si feedback | f_Si = 1 − w + w·h(Si)/h(Si_ref), h = Si/(Si + 5e-4 M) | Si_ref 2.7e-4 M, w 0.5 |
| `siup` rice Si uptake | − Si | 3.7 mol Si m-2 yr-1 · Monod([Si], 3e-4 M) |
| `nh4up` | − NH3 | §4.2 |

The Si feedback stands for silicon-strengthened roots and aerenchyma releasing more O2 (Fleck et al.
2011). It is the term that makes basalt lower CH4 in agreement with the field trials.

---

## 6. Gas transport and boundary schedules

These are boundary conditions and transport, kept separate from the chemistry so they can be matched
to ARTEMIS's own hydrology.

### 6.1 Gas transport (per cell, in-column sinks)

| Process | Form (mol m-3 yr-1) | Value | As a PHREEQC rate |
|---|---|---|---|
| `ch4pl` plant transport of CH4 | k(z)·φ·S_w·1000·K_H·(pCH4 − 1.9e-6), k(z) = ch4pl·exp(−z/0.1) | ch4pl 1100 yr-1 | first order on dissolved CH4 excess: k(z)/3.1536e7 · ([CH4] − K_H·1.9e-6) |
| `ch4rx` rhizosphere oxidation | Pox × the ch4pl rate is converted to CO2 in the cell; (1 − Pox) leaves to the atmosphere | Pox = pmin + (pmax − pmin)/(1 + ([H+]/H50)^n), `--pox=0.3:0.9:6:1` | |
| `n2opl` plant transport of N2O | as ch4pl, k scaled by n2opl | n2opl 0.207; atm N2O 330 ppb | |
| `ch4eb` ebullition | k·φ·S_w·1000·K_H·max(pCH4 − 0.4 atm, 0), vented directly | k 1000 yr-1 | |
| Diffusion | SCEPTER default: gas diffusion in air-filled pores plus aqueous diffusion, Millington-type tortuosity | | |

CH4 emission = surface diffusive flux + ebullition + (1 − Pox)·plant flux. The plant terms are zero in
the fallow.

### 6.2 Top boundary: floodwater exchange

While the surface is flooded, dissolved gases (O2, CO2, CH4, N2O, NH3) exchange between the top cell and
the atmosphere through the standing water with a transfer velocity:

Flux (mol m-2 yr-1) = k_fw·w_fw·1000·([X]_top − K_H·p_atm), k_fw = 36.5 m yr-1 (0.1 m d-1),
w_fw = min(1, max(0, (S_w,top − 0.9)/0.1)).

It replaces the half-cell diffusion at the top. When the surface drains (S_w < 0.9) gas-phase
diffusion through the air-filled pores takes over. Solutes enter with irrigation water at the
percolation rate (§6.4).

### 6.3 Flood, drain and fallow schedule

Water saturation in each cell is set by a surface saturation S_top(t) and a water-table depth z_wt:
S_w(z) = min(1, S_top + (1 − S_top)·z/z_wt). The surface saturation follows a yearly time series
(`Wet_temp.in`, linear interpolation, repeated every year).

| Phase | Timing | S_top | z_wt |
|---|---|---|---|
| Continuously flooded season | day 0 to 99 | 1 | — (all saturated) |
| Pre-harvest drain (`--phd=14`) | day 99 to 113 | 0.6 | 0.15 m |
| AWD season (`--awd`) | flooded day 0–29; then repeated cycles of drain 7.3 d (S_top 0.6, z_wt 0.15 m) and reflood 14.6 d, until the pre-harvest drain; ramps take 1.1 d | 1 / 0.6 | 0.15 m |
| Fallow (crop-free) | day 113 to 365 (0.69 yr) | 0.5 | 0.5 m |

Irrigated archetypes use the flooded schedule; rainfed archetypes use AWD. Tillage mixes the top 0.2 m
before each season (`TILL=0.2`, `paddy/multiyear/run_chain.sh`).

### 6.4 Percolation and irrigation water

| Item | Setting |
|---|---|
| Downward advection q | 0.5 m yr-1 (1.4 mm d-1), constant through the year |
| Column | 0.3 m, porosity 0.5, 15 cells in production (driver default 30) |
| Bottom | free outflow (advective leaching) |
| Irrigation / infiltrating water | Ca 0.5, Mg 0.2, Na 0.3, K 0.05, Si 0.2, SO4 0.2 mM, NO3 1e-3 mM, Cl 1.15 mM (+ charge balance) |

---

## 7. pH: what SCEPTER leaves to charge balance

SCEPTER solves pH in every Newton iteration from the electroneutrality of all total concentrations
(subroutine `calc_pH`), with speciation of DIC, total ammonia, sulfide, sulfate, Fe(II), Fe(III), Al,
Si and the exchanger. No reaction carries H+ in its stoichiometry. Each reaction changes only the
totals in §4 (for example `g2fe3` adds 4 Fe(II) and 1 DIC and removes 4 ferrihydrite). The 8 H+
consumed per CH2O is not written anywhere: charge balance raises pH to accommodate the added Fe2+.

So in PHREEQC, write the balanced forms of §4 (or give each KINETICS reaction as a neutral formula
of elements) and let PHREEQC compute pH. pH should then come out the same, provided the redox states
are decoupled as described in §2. The only places pH enters the rate laws are explicit factors:

| Reaction | pH factor |
|---|---|
| `fe2o2` | ([OH-])² (Singer–Stumm) |
| `fe2gr` | Ω contains [H+]² |
| `fe2fs` | Ω contains [H+] |
| `amo2o`, `g2n21`, `g2n0` | exp(−0.5(pH − 7)²) |
| N2O yield r | logistic, pH50 6.2 |
| `g2n2p` | 1/(1 + ([H+]/10^−5.5)²) |
| Pox (`ch4rx`) | logistic, pH50 6 |
| optional `--ch4ph` (methanogen bell curve), `--phdec`, `--phmo` | off in production; the ch4ph case is rejected by field data |

---

## 8. Temperature dependence

| Process | Treatment |
|---|---|
| `g2ch4` | Q10 = 3 |
| `g2fe3`, `g2gr`, `g2gt`, `g2so4`, `ch4o2`, `g2o2`, `amo2o`, `g2n21`, `g2n0`, `g2n2p`, `rtco2` | Q10 = 2 |
| Reference temperature | 25 °C: k(T) = k25·Q10^((T − 25)/10) (`paddy_q10`; overridable per reaction via `q10.in`) |
| Abiotic reactions (`fe2o2`, `fe2gr`, `groxy`, `fe2fs`, `fesox`, `sidox`), gas transport, `rtoxy`, `nh4up`, `siup` | no explicit temperature term |
| Mineral dissolution (basalt, calcite) | SCEPTER standard Arrhenius, activation energies from Palandri & Kharaka (2004) |
| Equilibrium constants, Henry's constants | SCEPTER thermodynamic temperature dependence |

Soil temperature is constant within a run (`--temp`, default 25 °C; archetype runs use the local
growing-season mean). The Q10 values are literature values (Segers 1998 range), not fitted.

---

## 9. Calibration: targets and data per rate

All targets are read on the **settled cycle**: a 3-year drained spin-up, then season + fallow years
(`run_chain.sh`, OM 450 g C m-2 yr-1 with 70% in season) until year 3–4, at the reference site
(25 °C, reference soil §10, flooded with pre-harvest drain, 120 kg urea-N). Re-check them in that order
after porting, because OM and Fe pools take years to settle.

| Constant | Value | Target and data | Model result (v15) |
|---|---|---|---|
| `g2ch4` vmax | 1630 yr-1 | IPCC 2019 flooded-rice baseline emission factor 1.19 kg CH4 ha-1 d-1 (≈ 134 kg per 113-d season) | 132–135 kg CH4 ha-1 season-1 |
| `g2ch4` Fe(III) Ki | 30 mol m-3 | CH4 onset after ~2–4 weeks of flooding | onset 10–25 d |
| `ch4pl` | 1100 yr-1 | plant-mediated share of emission ~80–90%, ebullition ~10–20% (field) | ebullition ~0% (reference) to 23% (acid high-OM archetype) |
| Pox (pmin, pmax, pH50, n) | 0.3, 0.9, 6, 1 | rhizosphere oxidation 0.4–0.9 in the literature; tested against steel-slag (Wang et al. 2014) and liming trials | root-zone pH 5.5 → Pox 0.47 |
| `ch4eb` | 1000 yr-1 | no post-harvest venting pulse; fallow CH4 ~10% of season | ~11% |
| `n2opl` | 0.207 × ch4pl | keeps absolute plant N2O conductance of v13 | |
| `g2n2p` vmax, pH50 | 300 yr-1, 5.5 | season N2O 0.3–0.4 kg N2O-N ha-1 with 120 kg N (EF ~0.3%, Akiyama et al. 2005); acid soils not > ~2.3 kg | 0.36 kg N2O-N |
| N2O yield r (n2oy) | 6.2, 1, 0.1 | N2O/(N2O+N2) versus pH from the literature (§4.2) | |
| O2 inhibition Ki | 2.5e-7 M | fallow N2O ~1 kg N ha-1 yr-1 (2.5e-5 gave ~7) | 0.3–0.9 kg N |
| `amo2o`, `g2n21` | 30, 300 | Maggi et al. 2008 defaults; not tuned | |
| `g2fe3`, `g2gt` | 300, 10 | not tuned individually; g2gt from crystalline/amorphous Fe reduction ratio | |
| f_gr (green-rust bioavailability) | 0.3 | pH 6.5–7 all season; porewater Fe2+ 1–3 mM (Ponnamperuma 1972) | pH 6.5–6.9 |
| `rtco2` | 60 mol C m-2 yr-1 | porewater pCO2 0.1–0.3 atm, HCO3- ~10–15 mM (Kirk et al. 2019) | |
| `rtoxy` | 3.65 mol O2 m-2 yr-1 | root radial O2 loss ~10 mmol m-2 d-1 (literature) | |
| Si feedback w | 0.5 | Gong et al. 2026: 36 t ha-1 basalt, AWD, CH4 −18% | −20% (best fit w ≈ 0.4) |
| `siup` | 3.7 mol m-2 yr-1 | rice Si uptake ~100–150 kg Si ha-1 season-1 | ~130 kg |
| Rock rate scaling | mafic ×0.45, plagioclase ×0.25, amorphous Si ×30 | Uchibayashi et al. 2025 pot fluxes of Mg and Na | match |
| `fe3in` | 0.3 mol Fe m-2 yr-1 | control-plot CH4 drift < ~1–2% per year over 20 yr | ~1.4%/yr |
| q | 0.5 m yr-1 | typical clay paddy percolation 1–2 mm d-1 | |
| AWD (not tuned, check) | | IPCC AWD scaling 0.55–0.77 | 0.66–0.69 |

### Independent validation (not used for tuning)

| Study | Treatment | Observed CH4 change | Model v15 |
|---|---|---|---|
| Gong et al. 2026 | basalt 36 t ha-1, AWD | −18% | −20% (w 0.5) |
| Park et al. 2026 | anoxic batch, basalt 10/20 wt% | −25 / −37% | −15 / −28% |
| Wang et al. 2014 | steel slag 2/4/8 t ha-1 | −27 / −43 / −49% | −24 / −27 / −30% |
| Wang F. et al. 2024 | wollastonite 5 t ha-1, soil inorganic C | +1.2 t CO2 ha-1 | +1.3 (calcite seeded) |
| Five liming trials | 1–4 t ha-1 lime | −20 to −53% | year 1 −8%, year 2 +5% (**misfit**) |

---

## 10. Soil and inputs (reference site)

| Input | Setting |
|---|---|
| OM | 450 g C m-2 yr-1 into the top 0.15 m; 70% in season, 30% in fallow |
| Fe(III) resupply | 0.3 mol Fe m-2 yr-1 as ferrihydrite with the OM input |
| Soil Fe(III) | ferrihydrite and goethite pools from the site (global proxy 0.012·clay/30 wt%, clipped 0.004–0.03) |
| CEC | 15 cmol kg-1; SCEPTER exchange log K: Na 5.9, K 4.8, Ca 10.47, Mg 10.786, Al 16.47; Fe2+ = Mg + 0.158 |
| Amorphous Si | 1% phytolith silica |
| Basalt | applied in the first 0.01 yr, mixed into the top 0.15 m, particle radius 10 µm |

---

## 11. Notes for re-implementation in PHREEQC

* Each microbial reaction maps to a `RATES` block with the Monod and inhibition product above; OM and
  Fe(III) minerals as `KINETICS` reactants. Evaluate the rate in SCEPTER units and convert (§1.2).
* Keep redox decoupled (§2), otherwise PHREEQC's pe equilibrium overrides the kinetic redox ladder.
* `fe2gr` and `fe2fs` are near-equilibrium pH-stats; they can be equilibrium phases once green rust has
  a PHASES entry, but the kinetic form avoids instantaneous Fe(II) uptake during drains.
* Gas transport (§6.1) has no PHREEQC equivalent in `TRANSPORT`; implement it as first-order removal of
  dissolved CH4 and N2O in each cell, and the floodwater exchange as a top-cell exchange term.
* Rate constants were calibrated inside SCEPTER's transport and grid. Re-run the §9 targets after porting.
