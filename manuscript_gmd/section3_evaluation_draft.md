## 3 Model evaluation against amendment experiments

### 3.1 Approach

We evaluated Riceland SCEPTER against five published experiments that added silicate rock or slag to paddy soil (Table 1). Together they test the model's main responses: CH₄ production, oxidation and emission (Wang et al., 2014; Gong et al., 2026; Park et al., 2026), mineral dissolution and cation release (Uchibayashi et al., 2025), and carbonate formation (Wang et al., 2024).

All runs use one frozen executable and one parameter set: the calibrated constants of the reference configuration, with the Si root feedback at K_Si = 2 mM and w = 0.2. Nothing was re-tuned for an individual site. Each site was run in three steps:

1. A 3-year drained spin-up of its own soil.
2. A 7-year settling period of season-plus-fallow cycles.
3. One experimental season, with control and treatment started from the same settled state.

Field sites use the 0.5 m column with a plough pan. Pot and batch experiments use a 0.3 m column without a pan, with the amendment mixed through the whole column. Soil acidity was set from the reported dry-soil pH:

- Calcareous soils (pH ≥ 6.8) received 0.5 wt% soil calcite.
- Other soils were given an irrigation-water alkalinity interpolated from pH: 0, 0.2, 0.8 and 2 meq L⁻¹ at pH 5.0, 6.0, 6.2 and 6.5.

Reported soil pH values are measured on dried slurries, not flooded porewater. The settled flooded porewater sits at pH 6.0–6.55, and calcite equilibrium with flooded-soil CO₂ caps it at about 6.55.

Two of the five tests are not fully independent of the calibration:

- The Si-feedback constants were chosen to minimise the joint misfit to Gong et al. (2026) and the 2 t ha⁻¹ treatment of Wang et al. (2014). The comparison for the 4 and 8 t ha⁻¹ treatments, and for CH₄ production and oxidation, is therefore a prediction, but the emission response at 2 t ha⁻¹ is not.
- The rock dissolution scaling factors (mafic minerals ×0.45, plagioclase ×0.25) were set against an earlier version of the Uchibayashi et al. (2025) release rates. That case is a consistency check on the current model, not an independent test.

We therefore present Section 3 as model evaluation, not independent validation. Table 3 marks each comparison as calibration-informed or as a prediction of a measurement not used in fitting. We keep the unresolved discrepancies (the ferric-Fe rise in Wang et al., 2014, and the carbonate gain in Wang et al., 2024) as identified targets for model development. We do not require the current version to reproduce every observation.

**Table 1.** Five experimental tests used to evaluate Riceland SCEPTER.

| Study | System and site | Amendment and application rates | Duration | Measured endpoints used | Model processes tested |
|---|---|---|---|---|---|
| Wang et al. (2014), Wetl. Ecol. Manag. | Field plots, subtropical China; continuous flooding (5 cm), final drain 2 weeks before harvest | Steel slag (SiO₂ 40.7, CaO 34.9 wt%), granular 1–7 mm; 0, 2, 4, 8 t ha⁻¹, applied 1 week before transplanting | One rice season, April–July 2011 (113 d modelled) | CH₄ emission (closed chambers); CH₄ production and oxidation (in-situ C₂H₂ inhibition, fortnightly); soil ferric Fe, available SiO₂ and pH at harvest | Si root feedback; partitioning of CH₄ production and oxidation; Fe(III) reduction and re-oxidation |
| Gong et al. (2026), J. Environ. Manag. | Planted pots (30 cm diameter × 50 cm, 35 kg soil), calcareous soil (pH ~7.4); flooding–drainage–reflooding–moist cycle; 240 kg N ha⁻¹; n = 3 | Basalt 2% w/w (~36 t ha⁻¹) mixed into the upper 10 kg of soil (biochar treatments not modelled) | 132 d to harvest (131 d modelled) | Seasonal CH₄ and N₂O emission; soil carbonate (inorganic C) | Si root feedback with plants and intermittent drainage; calcareous-soil carbonate system |
| Park et al. (2026), Appl. Biol. Chem. | Anoxic batch incubations of Korean paddy soil, no plants, 25 °C, dark; biotic and autoclaved, n = 3 | JB-1b basalt at 0, 10, 20 wt% of soil (390 and 780 t ha⁻¹ equivalent in the 0.3 m model column) | 100 d for gases; cations followed for 2 years | Headspace CH₄ and CO₂; dissolved Si, Ca, Mg; DOC; SOC; methanogen community | Basalt geochemistry (Fe(III) supply, pH, cations) without plant transport, root O₂ or the Si feedback |
| Uchibayashi et al. (2025), Soil Sci. Plant Nutr. | Planted pots, 8 kg Gray Lowland soil (pH 6.19); continuous flooding | Akita basalt (bytownite 37, labradorite 26, pyroxene 10, magnetite 7, smectite 6, olivine 4 wt%; D50 120 µm) at 0, 5, 10, 20, 50, 100 t ha⁻¹ | 99 d | Release of Si, Ca, Mg, K and Na by mass balance; exchangeable cations; available Si; pH; soil C. No gas fluxes | Mineral dissolution kinetics under flooding; cation release by mineral |
| Wang et al. (2024), Plant Soil | Field, Liaoning, NE China; eight 14 × 11 m plots, n = 4; intermittent flooding | Wollastonite 5 t ha⁻¹, tilled to 15 cm after irrigation, before transplanting | One season, May–October 2022 (153 d modelled) | Soil inorganic C at 0–20 and 20–40 cm (every 2 months); pH; available Si; yield and grain quality; grain As and Cd | Carbonate precipitation; Ca partitioning among calcite, the exchange complex and leaching; percolation |

Results are summarised in Table 3. Each case study below covers four things in turn: the experimental conditions, the model configuration, observed against simulated values, and the mechanistic interpretation and its limits.

### 3.2 Wang et al. (2014): CH₄ suppression by steel slag

**Experimental conditions.**
- Steel slag was applied at 2, 4 and 8 t ha⁻¹ to continuously flooded field plots in subtropical China.
- Seasonal mean CH₄ emission fell from 3.11 to 2.29, 1.76 and 1.59 mg m⁻² h⁻¹ (−27, −43, −49%).
- CH₄ production, measured by in-situ C₂H₂ inhibition, fell by 16, 33 and 41%.
- Oxidation stayed nearly constant (0.46–0.59 mg m⁻² h⁻¹).
- At harvest, ferric Fe rose from 3.2 to 8.8 g kg⁻¹ at 8 t ha⁻¹, available SiO₂ rose from 254 to 1232 mg kg⁻¹, and soil pH rose from 6.48 to 7.16.
- The authors attributed the CH₄ decrease to the slag's Fe oxides acting as electron acceptors, and to higher soil redox.

**Model configuration.**
- Field column with plough pan, continuous flooding, final drain 14 d before harvest, 113 d season.
- Non-calcareous soil (pH 6.5) with 2 meq L⁻¹ irrigation alkalinity.
- Slag represented by its oxide composition, with a lognormal grain-size distribution matching the reported 1–7 mm granules (volume median 2.65 mm).
- Model production is methanogenesis. Model oxidation is bulk-soil oxidation plus rhizosphere oxidation, which is the flux that in-situ C₂H₂ blocks.

**Observed versus simulated.**
- The simulated control emitted 2.95 mg CH₄ m⁻² h⁻¹ (observed 3.11).
- Emission fell by 35, 43 and 52% (observed 27, 43 and 49%).
- Production fell by 27, 33 and 38% (observed 16, 33 and 41%).
- Simulated oxidation fell with production, by 24, 29 and 33%, so the oxidised fraction stayed at 73–79%. It did not rise.
- Reactive Fe(III) at harvest in 0–15 cm rose by only 0.5 g kg⁻¹ at 8 t ha⁻¹, against an observed rise of 5.6 g kg⁻¹ in extractable ferric Fe.
- In the earlier model version (v21, otherwise identical physics), switching off the Si feedback (w = 0) left slag changing emission by less than 10% at every dose.

**Mechanistic interpretation and limitations.**

The model reproduces the size and dose response of the CH₄ decrease, and it does so through lower production rather than more oxidation, as observed. We regard this as a successful prediction of the CH₄ response. It is not validation of the mechanism. In the model the decrease runs through the Si root feedback:

- Higher porewater Si raises radial O₂ loss from roots (1.24 to 2.29 mol O₂ m⁻² per season at 8 t ha⁻¹).
- The extra O₂ is consumed mainly by aerobic respiration of organic matter.
- That respiration takes substrate away from both Fe(III) reduction (−32%) and methanogenesis.

Fe(III) reduction and Fe(II) re-oxidation therefore both *fall* in the treatment. That is the opposite of the pathway the authors proposed: more root O₂, more Fe(II) oxidation, more reactive Fe(III), less methanogenesis.

We tested that pathway explicitly by raising the catalysed Fe(II) oxygenation rate by six to seven orders of magnitude while conserving the root O₂ supply. This was a labelled sensitivity, not a change to the reference. It failed all four diagnostic tests (Fe(III) gain, Fe(III) reduction, production and oxidation). The reason is that the root O₂ increase (+0.75 mol O₂ m⁻²) could oxidise at most 3 mol Fe m⁻², about 1 g kg⁻¹ over 15 cm, an order of magnitude less than the observed ferric-Fe gain.

The observed ferric-Fe rise therefore remains unexplained. It cannot come from the slag itself. The slag's Fe (4.8 wt% Fe₂O₃; 33.6 g Fe kg⁻¹, or 269 kg Fe ha⁻¹ at 8 t ha⁻¹, as reported by Wang et al.) is included in the model as reactive ferrihydrite. Spread through the 0–30 cm sampled layer (bulk density ~1.0 g cm⁻³, ~3000 t soil ha⁻¹), it adds only ~0.09 g Fe kg⁻¹. The observed increases are about 50 times larger: +5.6 g kg⁻¹ ferric Fe and +4.4 g kg⁻¹ 1 M HCl-extractable total Fe, equivalent to ~13–17 t Fe ha⁻¹. They therefore reflect a change in how much native soil Fe the extraction recovers, or sampling variability, not an Fe input. Candidate explanations:

- Changes in Fe(II) speciation or extractability at higher pH and Ca.
- Mn oxides.
- The extraction method. Extractable "ferric Fe" is not equivalent to the model's microbially reducible Fe(III) pool.

We report this as a model limitation. We have not tuned it away.

A second limitation is the oxidised fraction. The model oxidises about 75% of produced CH₄, against 10–35% implied by Wang's C₂H₂ data. Field C₂H₂ inhibition is known to underestimate oxidation, but the gap is large.

### 3.3 Gong et al. (2026): basalt in planted pots

**Experimental conditions.**
- Basalt at about 36 t ha⁻¹ was mixed into the upper part of planted pots of a calcareous paddy soil.
- Water was managed through flooding, mid-season drainage, reflooding and a moist phase.
- Seasonal CH₄ fell from 52.4 ± 6.5 to 43.1 ± 3.1 kg ha⁻¹ (−18%, n = 3, not significant at p < 0.05).
- The authors also reported large soil carbonate formation, about 11.3 t CO₂ ha⁻¹.

**Model configuration.**
- 0.3 m pot column without a pan. Basalt is mixed into the top 0.086 m, the same fraction of the column (10/35) as the upper 10 kg of the 35 kg pot.
- Calcareous soil (0.5 wt% calcite).
- Intermittent drainage (AWD) and 240 kg N ha⁻¹.
- 131 d season.
- Basalt mineralogy from the paper's XRD, with an 80–100 µm grain size.

**Observed versus simulated.**
- The simulated control emitted 85 kg CH₄ ha⁻¹ (observed 52), about 60% too high.
- Basalt lowered CH₄ by 23% (85 to 65 kg ha⁻¹; observed −18%). Mixing the same dose through the whole column instead gave −25%, so amendment placement changes the result by only 2 points.
- With the Si feedback switched off, the change was −10%.
- Simulated calcite gain was 1.1 t CO₂ ha⁻¹, from rock cation release equivalent to 3.5 t CO₂ ha⁻¹.

**Mechanistic interpretation and limitations.**

The model overpredicts the relative CH₄ decrease by about 30% (−23% against −18%), and the control emission by about 60%. We state this explicitly, because together they make Gong et al. the weakest CH₄ comparison in the evaluation set. Gong et al. and the 2 t ha⁻¹ slag treatment of Wang et al. (2014) reach similar simulated porewater Si (~1.3 mM), yet the observed CH₄ cut is about twice as large for the slag. No single Si-response curve fits both. The adopted constants (K_Si 2 mM, w 0.2) minimise the joint misfit.

Possible reasons for the gap:

- The Gong pots were drained repeatedly. Drainage already aerates the root zone, so the extra benefit of root O₂ loss may be smaller than under continuous flooding.
- The calcareous soil buffers the pH rise that the basalt would otherwise cause.
- The observed effect is itself uncertain (n = 3, p ≈ 0.1).

Without the Si feedback, basalt geochemistry alone lowers CH₄ by 10%, about half of the observed effect.

Gong et al. report a carbonate gain of 11.3 ± 3.1 t CO₂ ha⁻¹. That is three times the CO₂ equivalent of all the cations the model releases (3.5 t CO₂ ha⁻¹). It is also larger than the total Ca + Mg in 36 t ha⁻¹ of a typical basalt (~5 t CO₂ ha⁻¹ at 9 wt% CaO and 7 wt% MgO). We therefore do not use it as a quantitative target.

### 3.4 Park et al. (2026): basalt without plants

**Experimental conditions.**
- Anoxic, dark, 25 °C batch incubations of a Korean paddy soil with 10 and 20 wt% JB-1b basalt.
- CH₄ appeared after 50–60 days in the biotic batches.
- By day 100, cumulative CH₄ was about 75% and 63% of the control.
- Autoclaved batches produced no CH₄ but released similar Si, Ca and Mg.

**Model configuration.**
- 0.3 m closed column with negligible percolation (q = 10⁻⁴ m yr⁻¹).
- No plants: plant transport, root respiration, root O₂ loss, N uptake and Si uptake are all off, so the Si root feedback cannot act.
- 2 meq L⁻¹ water alkalinity.
- JB-1b mineralogy, with 3% of rock Fe as Fe(III).
- Basalt mixed through the column at 390 and 780 t ha⁻¹ equivalent.
- 100 d.

**Observed versus simulated.**
- CH₄ fell by 20% and 34% at 10 and 20 wt% (observed about −25% and −37%).

**Mechanistic interpretation and limitations.**

Park et al. is the one test in which the Si root feedback is absent by construction. It isolates the soil geochemistry: rock Fe(III) acting as an electron acceptor, the higher pH, and cation supply. The model reproduces both the size and the dose scaling of the decrease. In the model's electron budget, basalt raises Fe(III) reduction by 0.54 and 0.94 mol C m⁻² at 10 and 20 wt%, while methanogenesis falls by 1.57 and 2.67 mol C m⁻². Electron-acceptor competition from rock Fe(III) therefore accounts directly for only about a third of the decrease. The rest is lower total organic-matter mineralisation, which in the model arises because Fe(III) inhibits methanogenesis without fully replacing it as an electron sink. Both terms depend on the rock's reactive Fe(III), which was assumed (3% of rock Fe), not measured. We therefore treat the mechanism as uncertain and do not claim that basalt-derived Fe(III) is the established cause of the observed suppression.

At field doses, the same geochemistry gives a small effect (Sections 3.2 and 3.3). The large field responses therefore need the plant pathway. The test is limited in three ways:

- The doses are 10–40 times typical field rates.
- The decrease depends on the assumed rock Fe(III) fraction (3%), which was not measured.
- The batch CH₄ lag (50–60 d) reflects microbial community development that the model does not represent.

### 3.5 Uchibayashi et al. (2025): dissolution and cation release

**Experimental conditions.**
- Akita basalt at 0–100 t ha⁻¹ was mixed into continuously flooded, planted pots of a Gray Lowland soil (pH 6.19) for 99 days.
- Cation release was computed by mass balance.
- At 100 t ha⁻¹, Mg was released at 6.5 and Na at 4.6 mmol m⁻² d⁻¹.
- Ca release did not scale with dose.
- Soil pH (6.5–6.7 after harvest) and soil C did not change significantly with dose.
- Gases were not measured.

**Model configuration.**
- 0.3 m pot column, continuous flooding, 99 d.
- Irrigation alkalinity 0.8 meq L⁻¹, matching pH 6.2.
- Akita mineralogy and the measured grain-size distribution (D10/D50/D90 = 22/120/370 µm).
- Basalt at 100 t ha⁻¹ mixed through the column.

**Observed versus simulated.**
- Simulated Na release was 4.8 mmol m⁻² d⁻¹ (observed 4.6).
- Simulated Mg release was 12.1 mmol m⁻² d⁻¹ (observed 6.5).

**Mechanistic interpretation and limitations.**

Na comes almost entirely from plagioclase. Its agreement shows that the plagioclase rate, scaled ×0.25 from laboratory rates, is consistent with these pots. Mg is released about twice too fast even with the mafic minerals scaled ×0.45.

The mineral-by-mineral dissolution fluxes locate the excess. Olivine supplies 11.4 of the 12.1 mmol Mg m⁻² d⁻¹ (94%), and augite the remaining 0.65. About 19% of the olivine dissolves in 99 days, against 16% of the bytownite and under 2% of the labradorite and augite. Olivine makes up only 4 wt% of the rock, but its rate law is two to three orders of magnitude faster than pyroxene's. With dissolution close to linear at this extent, an olivine rate about half the current value (≈0.5×, applied to olivine alone) would bring Mg release to the observed 6.5 mmol m⁻² d⁻¹ without changing Na. The olivine explanation is therefore quantitatively sufficient. We do not adopt that adjustment in the reference parameterisation. Possible explanations:

- The rock's olivine may be partly altered (iddingsite, serpentine) and less reactive than fresh olivine.
- Released Mg may be retained in secondary smectite or on exchange sites, and so be missed by the mass balance.
- The model's single scaling factor for all mafic minerals cannot treat olivine and pyroxene separately.

Because the scaling factors were set against these data, the case shows internal consistency rather than independent skill. It also shows that a single mafic factor is too coarse for olivine-bearing rocks.

### 3.6 Wang et al. (2024): carbonate accumulation under wollastonite

**Experimental conditions.**
- Wollastonite at 5 t ha⁻¹ was tilled into 0–15 cm of intermittently flooded field plots in Liaoning.
- Over one season, soil inorganic C in 0–20 cm rose by 1.20 ± 0.03 t CO₂ ha⁻¹ relative to the control.
- pH rose from 6.85 to 7.24.
- The authors interpreted the inorganic-C gain as CO₂ sequestered as pedogenic carbonate.

**Model configuration.**
- Field column with plough pan.
- Calcareous soil (0.5 wt% calcite).
- Intermittent drainage with a pre-harvest drain, 153 d season.
- Wollastonite (with its 2.6 wt% CaCO₃) tilled into 0–15 cm, about 3.7 mol Ca m⁻².
- No percolation data are reported. We ran three constant rates, 1.1, 0.5 and 0.25 m yr⁻¹, applied uniformly through the season.

**Observed versus simulated.**

| Treatment minus control, 0–20 cm (mol m⁻² unless stated) | q 1.1 m yr⁻¹ | q 0.5 | q 0.25 |
|---|---|---|---|
| Calcite gain (t CO₂ ha⁻¹) | 0.31 | 0.40 | 0.48 |
| Ca exported through 20 cm | 0.60 | 0.24 | 0.05 |
| Exchangeable Ca | +2.35 | +2.49 | +2.46 |
| Season-mean log₁₀ Ω_calcite, control → treatment | −1.18 → −0.35 | −1.00 → −0.20 | −0.90 → −0.13 |
| Season CH₄ change | −46% | −48% | −48% |

**Calcium conservation.** The whole-column Ca budget, taken from the model's own time-integrated fluxes, closes to better than 10⁻¹¹ mol m⁻² in every run. In treatment minus control (mol m⁻²):

| Term | q 1.1 m yr⁻¹ | q 0.5 | q 0.25 |
|---|---|---|---|
| Ca added: wollastonite + its native calcite | 3.587 + 0.131 = 3.718 | 3.718 | 3.718 |
| Wollastonite remaining | 0.002 | 0.002 | 0.002 |
| Exchangeable Ca | +2.77 | +2.75 | +2.59 |
| Dissolved Ca (porewater) | +0.10 | +0.10 | +0.09 |
| Calcite, whole column | +0.90 (0.39 t CO₂ ha⁻¹) | +0.91 (0.40) | +1.06 (0.47) |
| Ca leached out of the column base | −0.05 | −0.05 | −0.02 |
| Sum of sinks | 3.72 | 3.72 | 3.72 |

There are no other Ca sinks in the model: no plant Ca uptake and no other Ca-bearing secondary phase. The Ca that passes below 20 cm (0.05–0.60 mol m⁻², see the percolation table below) is retained deeper in the column, mainly on exchange sites. Leaching from the column base is slightly lower in the treatment than in the control. The carbonate shortfall is therefore not an accounting error. In the model, the added Ca goes mainly to the exchange complex rather than to calcite.

**Mechanistic interpretation and limitations.**

The wollastonite dissolves completely in all three cases. About two thirds of its Ca ends up on the exchange complex. Calcite forms only in short supersaturated windows during drainage, because the plough layer is undersaturated on season average. Lower percolation shifts Ca from leaching into calcite. Even at q = 0.25 m yr⁻¹, however, simulated carbonate reaches only 0.48 t CO₂ ha⁻¹, 40% of the observed 1.20.

We did not fit percolation to the observation. Wang et al. give no hydrological data that would independently justify a low rate. The seasonal pattern of percolation (ponded seepage during flooding, none during drains) is also not resolved by the model's constant rate, and this may matter for calcite, which forms during drains.

The observed gain corresponds to about 91% of the wollastonite Ca precipitating as carbonate within one season. That is hard to reconcile with the exchange capacity of the soil. It may also include carbonate redistributed from the control or measurement variability, since it is a difference between two plot sets.

The observed treatment–control difference in carbonate stock should not be equated with atmospheric CO₂ removal, for three reasons:

- Part of the carbonate C may come from soil respiration rather than the atmosphere.
- Carbonate formed in the plough layer can redissolve in later seasons. In multi-year runs, calcite formed during the season partly redissolves during the fallow (Section 4).
- Ca leached as bicarbonate is CDR only once it reaches drainage water without degassing.

We therefore treat carbonate stock as an observable and report CDR separately, as potentially durable alkalinity exported below the root zone (Section 4).

**Table 3.** Summary of observed and simulated responses (frozen model, K_Si 2 mM, w 0.2).

| Test | Metric | Observed | Model | Assessment | Role in calibration |
|---|---|---|---|---|---|
| Wang 2014 | CH₄ emission change at 2 / 4 / 8 t ha⁻¹ | −27 / −43 / −49% | −35 / −43 / −52% | good (2 t ha⁻¹ used in Si calibration) | calibration-informed at 2 t ha⁻¹; prediction at 4, 8 t ha⁻¹ |
| Wang 2014 | CH₄ production change | −16 / −33 / −41% | −27 / −33 / −38% | good | prediction |
| Wang 2014 | CH₄ oxidation | ~constant | −24 / −29 / −33% | oxidation does not rise, as observed; absolute fraction too high | prediction |
| Wang 2014 | ferric Fe gain at 8 t ha⁻¹ | +5.6 g kg⁻¹ (extractable) | +0.5 g kg⁻¹ (reactive Fe(III)) | fails; mechanism not reproduced | prediction |
| Gong 2026 | CH₄ change at 36 t ha⁻¹ | −18% (p ≈ 0.1) | −23% (−10% without Si feedback) | overpredicted ~1.3× (used in Si calibration); control 85 vs 52 kg ha⁻¹ | calibration-informed |
| Park 2026 | CH₄ change at 10 / 20 wt% | ~−25 / −37% | −20 / −34% | good (no plant pathway) | prediction |
| Uchibayashi 2025 | Na / Mg release at 100 t ha⁻¹ (mmol m⁻² d⁻¹) | 4.6 / 6.5 | 4.8 / 12.1 | Na good; Mg ~2× too fast | calibration-informed |
| Wang 2024 | inorganic C gain, 0–20 cm (t CO₂ ha⁻¹) | 1.20 ± 0.03 | 0.31–0.48 (q 1.1–0.25 m yr⁻¹) | underpredicted ~2.5–4× | prediction |

**Table 4.** Reproducibility of the evaluation runs. All runs use executable `scepter_paddy_v22` (commit 162913a) with `krxn_calibrated` constants (methanogenesis vmax 1860 yr⁻¹) and the Si feedback `--rtsi=2e-3:2.7e-4:0.2`. Each run restarts from its site's settled state (3-yr drained spin-up + 7-yr settle).

| Test | Settled state | Runs (control / treatments) | Site options |
|---|---|---|---|
| Wang 2014 | field column with pan, irrigation alkalinity 2 meq L⁻¹ | z14_0 / z14_2, z14_4, z14_8 | 113 d, continuous flooding, 14-d pre-harvest drain, slag PSD `psdrain_wangslag.in` |
| Gong 2026 | 0.3 m column, 0.5 wt% calcite (vc_r03c05) | pgongc_0 / pgongd_36 (top 0.086 m); pgongc_36 (whole column) | 131 d, AWD, 240 kg N ha⁻¹ |
| Park 2026 | 0.3 m column, 2 meq L⁻¹ (vc_r03i10) | pparka_0 / pparka_390, pparka_780 | 100 d, no plant processes, q 10⁻⁴ m yr⁻¹, rock Fe(III) 3% |
| Uchibayashi 2025 | 0.3 m column, 0.8 meq L⁻¹ (vc_r03i4) | puchia_0 / puchia_100 | 99 d, continuous flooding, Akita PSD |
| Wang 2024 | field column with pan, 0.5 wt% calcite | z24_0 / z24_5 (q 1.1); z24q05_*, z24q025_* | 153 d, AWD, wollastonite tilled to 0.15 m |
