# Riceland SCEPTER: process specification (v15)

This document lists every process added to SCEPTER for flooded rice paddies, as rate laws with
stoichiometry, parameters, units, inhibition terms and the source of each value. It is written so
the processes can be re-implemented in another reactive-transport code (for example as PHREEQC
`RATES`/`KINETICS` blocks) without porting the Fortran. Code references are to `scepter.f90` and
`paddy/run_flooded_column.py` on branch `paddy/fe3-reduction` at the commit that adds this file
(model version v15, calibration `paddy/krxn_calibrated_v15.in`).

Conventions follow SCEPTER: time in years, depth in m, solid and gas concentrations in mol per m3
of bulk soil, aqueous concentrations in mol/L (M), partial pressures in atm. Reaction rates are in
mol per m3 bulk per yr unless stated. "OM" is the reactive organic-matter pool `g2` (CH2O
stoichiometry). "Monod(X, K)" means X/(X + K); "Inh(X, Ki)" means Ki/(X + Ki).

## 1. Common modifiers

| Modifier | Form | Value | Source |
|---|---|---|---|
| OM limitation | Monod(g2, 100 mol/m3) on all OM-consuming reactions | Km 100 mol/m3 | SCEPTER default |
| Saturation factor | f_sat = min(2·S_w, 1) on anaerobic reactions | | this work |
| O2 inhibition (anaerobic pathways) | Inh([O2], Ki) | Ki = 2.5e-7 M (`--o2inh`) | v7 calibration; the upstream 2.52e-5 left ~9% of anaerobic rates in aerated soil |
| NO3 inhibition (Fe, SO4, CH4 pathways) | Inh([NO3-], 11.3e-5 M) | | Maggi et al. 2008 |
| Temperature | Q10^((T − 25)/10) | Q10 = 3 methanogenesis, 2 other microbial reactions; Tref 25 °C (override via `q10.in`) | Segers 1998 range; this work |
| Plant factor | 1 in season, 0 in fallow (`plant_temp.in`) on all root/plant terms | | this work |

All microbial rates below are multiplied by the temperature factor and, where they consume OM, by
the OM Monod term.

## 2. Carbon, iron and sulfur reactions

| Reaction | Stoichiometry | Rate law | Parameters (v15) |
|---|---|---|---|
| `g2o2` aerobic respiration | CH2O + O2 → CO2 + H2O | v·Monod(g2)·Monod([O2], 5e-6 M) | v 300 /yr |
| `g2fe3` ferrihydrite reduction | CH2O + 4 Fe(OH)3 + 8 H+ → CO2 + 4 Fe2+ + 11 H2O | v·Monod(g2)·Fe_pool/(Fe_pool + K)·f_O2·f_NO3·f_sat, Fe_pool = amfe3 + 2·f_gr·GR | v 300 /yr, K 100 mol/m3 |
| `g2gr` green-rust reduction | CH2O + 2 GR → CO2 + 6 Fe2+ (charge balanced with OH-) | shares the `g2fe3` expression; GR enters Fe_pool with bioavailability f_gr | f_gr 0.3 |
| `g2gt` goethite reduction | as `g2fe3` with FeOOH | v·Monod(g2)·Monod(gt, 100 mol/m3)·f_O2·f_NO3·f_sat | v 10 /yr (crystalline Fe ~30x slower; Roden & Zachara 1996) |
| `g2so4` sulfate reduction | 2 CH2O + SO4 2- → 2 CO2 + HS- (+ H2O, H+ balance) | v·Monod(g2)·Monod([SO4], 1e-4 M)·f_O2·f_NO3·f_sat | v 300 /yr |
| `g2ch4` methanogenesis | 2 CH2O → CH4 + CO2 | v·Monod(g2)·Inh_Fe·Inh([SO4], 1e-4 M)·f_O2·f_NO3·f_sat, Inh_Fe = 30/(Fe_pool + 30) | v 1630 /yr, Fe(III) Ki 30 mol/m3, Q10 3 |
| `ch4o2` methanotrophy (bulk soil) | CH4 + 2 O2 → CO2 + 2 H2O | v·Monod([CH4], 5e-6 M)·Monod([O2], 2e-5 M) | v 500 mol/m3/yr |
| `fe2o2` Fe(II) oxidation | Fe2+ + 1/4 O2 + 5/2 H2O → Fe(OH)3 + 2 H+ | (k_SS·(Kw/[H+])² + k0)·[Fe2+]·pO2, applied to dissolved and exchangeable Fe(II) | k_SS 8e13 (Singer & Stumm 1970), k0 1e-7 /min/atm |
| `fe2gr` green-rust pH-stat | Fe2+ + 2 Fe(OH)3 → GR + 2 H+ | v·f_sat·ln Ω with Ω = 10^(−9.98)·[Fe2+]/[H+]²; forward × Monod(amfe3, 1), reverse × Monod(GR, 1) | v 1e4 |
| `groxy` green-rust oxidation | GR + 1/4 O2 → 3 Fe(OH)3 | 1e7·GR·k_H,O2·pO2 | |
| `fe2fs` FeS precipitation | Fe2+ + HS- → FeS + H+ | 1e3·ln Ω, log K 3.915 | |
| `fesox` FeS oxidation | FeS + 2 O2 → Fe2+ + SO4 2- | 1e6·FeS·[O2] | |
| `sidox` siderite oxidation | FeCO3 + 1/4 O2 + 3/2 H2O → Fe(OH)3 + CO2 | 1e5·siderite·[O2] | needed to stop siderite piling up below the plough layer |

Green rust stands in for the observed pH-buffering of flooded soils near pH 6.5–7 (Ponnamperuma
1972). Fe(II) sorption on ferrihydrite is a switchable alternative. Goethite and hematite are not
allowed to precipitate; calcite nucleation uses a seed of 0.001.

## 3. Nitrogen reactions

| Reaction | Stoichiometry | Rate law | Parameters (v15) |
|---|---|---|---|
| `amo2o` nitrification | NH4+ + 2 O2 → NO3- + 2 H+ + H2O | v·Monod([NH4], 14e-5 M)·Monod([O2], 2.41e-5 M)·exp(−0.5(pH − 7)²) | v 30 /yr; Maggi et al. 2008 |
| `g2n21` denitrification to N2O | CH2O + NO3- → N2O (stoichiometric fraction r) | v·Monod(g2)·Monod([NO3], 11.3e-5 M)·f_O2·exp(−0.5(pH − 7)²)·r | v 300 /yr |
| `g2n0` denitrification to N2 | as above to N2, rate 1.25·(1 − r)·(g2n21 rate) | | same vmax |
| N2O yield | r = 0.1 + 0.9/(1 + 10^(pH − 6.2)) (`--n2oy=6.2:1:0.1`) | | Simek & Cooper 2002; Liu et al. 2010; Russenes et al. 2016 |
| `g2n2p` N2O reduction | CH2O + 2 N2O → CO2 + 2 N2 | v·Monod(g2)·Monod([N2O], 5e-6 M)·Inh([NO3], 1e-3 M)·f_O2·1/(1 + ([H+]/H50)²) | v 300 /yr, pH50 5.5 |
| `nh4up` plant N uptake | NH4+ removed | 3.5 mol N/m2/yr·Monod([NH4], 5e-5 M), on root profile | |

Fertiliser: urea, 120 kg N/ha by default, spread evenly over the 113-day season. Residue N enters
with the OM rain at N:C 0.03 (`--n2c`).

## 4. Plant and root terms

All are distributed over the root profile (exponential, e-folding 0.1 m) and switched off in the fallow.

| Term | Form | Value |
|---|---|---|
| `rtco2` root respiration | CO2 source | 60 mol C/m2/yr |
| `rtoxy` radial O2 loss | O2 source × Si feedback f_Si | 3.65 mol O2/m2/yr (10 mmol/m2/d) |
| Si feedback | f_Si = 1 − w + w·h(Si)/h(Si_ref), h = Si/(Si + 5e-4 M) | Si_ref 2.7e-4 M, w 0.5 (Gong et al. fit ~0.4–0.53) |
| `siup` rice Si uptake | 3.7 mol Si/m2/yr·Monod([Si], 3e-4 M) | ~130 kg Si/ha/season without basalt |

The Si feedback represents silicon-strengthened rice roots and aerenchyma releasing more O2
(Fleck et al. 2011); it is the term that makes basalt lower CH4 in agreement with field trials.

## 5. Gas transport

| Process | Form | Value |
|---|---|---|
| `ch4pl` plant-mediated CH4 transport | k(z)·φ·S_w·1e3·K_H·(pCH4 − 1.9e-6 atm), k(z) = ch4pl·exp(−z/0.1 m) | ch4pl 1100 /yr |
| `ch4rx` rhizosphere oxidation of transported CH4 | a fraction Pox of the `ch4pl` flux is oxidised to CO2; the rest (1 − Pox) is emitted | Pox = pmin + (pmax − pmin)/(1 + ([H+]/H50)^n), `--pox=0.3:0.9:6:1` |
| `n2opl` plant N2O transport | as `ch4pl` with N2O, k scaled by 0.207 | |
| `ch4eb` ebullition | k·φ·S_w·1e3·K_H·max(pCH4 − 0.4 atm, 0), vented to the atmosphere | k 1000 /yr |
| Floodwater exchange | gas transfer at the soil surface with a standing water layer | k_fw 36.5 m/yr |
| Diffusion | SCEPTER default gas + aqueous diffusion with saturation-dependent tortuosity | |

Pox replaces the fixed 50% of the earlier versions; it lowers emission at near-neutral pH
(observed in lime and slag trials). Emission = surface diffusive flux + ebullition + (1 − Pox)·plant flux.

## 6. Water, flooding and the annual cycle

| Item | Setting |
|---|---|
| Column | 0.3 m, 15 cells, porosity 0.5 |
| Season | 113 days (0.31 yr), saturation 1 (flooded); 14-day pre-harvest drain (`--phd=14`) |
| AWD | drain to 0.15 m with surface saturation 0.6 on the Wet_temp schedule (`--awd`); rainfed sites use AWD |
| Fallow | 0.69 yr, surface saturation 0.5, water table 0.5 m, plant terms off |
| Percolation | q = 0.5 m/yr (1.4 mm/d) constant downward advection |
| Tillage | the top 0.2 m is mixed before each season (`TILL=0.2`) |
| Irrigation water | Ca 0.5, Mg 0.2, Na 0.3, K 0.05, Si 0.2, SO4 0.2 mM, Cl 1.15 mM plus charge balance |

## 7. Inputs and soil

| Input | Setting |
|---|---|
| OM | g2 rain into the top 0.15 m, 70% in season and 30% in fallow (`OMF=0.7`) |
| Fe(III) resupply | 0.3 mol Fe/m2/yr as ferrihydrite with the OM rain (`--fe3in`), balancing Fe(II) leaching |
| Reference soil | CEC 15 cmol/kg; exchange log K Na 5.9, K 4.8, Ca 10.47, Mg 10.786, Al 16.47; Fe2+ = Mg + 0.158; 1% phytolith amorphous Si (rate ×30) |
| Rock kinetics scaling | mafic minerals ×0.45, plagioclase ×0.25 (fit to Uchibayashi et al. column fluxes) |

## 8. Calibration summary (v15)

| Constant | Value | Target |
|---|---|---|
| g2ch4 vmax | 1630 /yr | IPCC 2019 flooded rice ~1.19 kg CH4/ha/d → 132–135 kg/season on the reference site |
| ch4pl | 1100 /yr | plant transport carries most emission; ebullition ~10–20% |
| n2opl | 0.207 | keeps absolute plant N2O conductance |
| g2n2p | 300 /yr, pH50 5.5 | ~0.36–0.4 kg N2O-N/season (EF ~0.3%, Akiyama et al. 2005) |
| ch4eb | 1000 /yr | stops a post-harvest venting pulse |
| f_gr | 0.3 | |
| Si feedback w | 0.5 (fit ~0.4) | Gong et al. basalt trial −18% CH4 |

Validation against field data (v15): Gong 36 t/ha basalt −20% (w 0.5) vs −18% observed; Park
silicate slag −15/−28% vs −25/−37%; Wang et al. 2014 steel slag −24/−27/−30% vs −27/−43/−49%;
lime trials not reproduced (model year 1 −8%, year 2 +5% vs −20 to −53% observed).

## 9. Notes for re-implementation in PHREEQC

* Each microbial reaction maps to a `RATES` block with a Monod/inhibition product; OM as a
  kinetic reactant `CH2O`, Fe(III) minerals as `KINETICS` reactants.
* `fe2gr`, `fe2fs` are near-equilibrium pH-stats and can be equilibrium phases.
* Gas transport (plant, ebullition, floodwater) needs explicit sink terms in each cell; PHREEQC
  `TRANSPORT` has no gas phase, so treat the three pathways as first-order removal of dissolved
  CH4/N2O.
* Rate constants were calibrated in this model's transport and should be re-checked after porting.
