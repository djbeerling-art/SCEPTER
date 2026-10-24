# Paddy (flooded rice soil) development

Work towards a riceland version of SCEPTER for enhanced weathering in rice paddies.

## Added so far

- `g2fe3` external reaction in `scepter.f90`: dissimilatory Fe(III) reduction coupled to OM oxidation,
  CH2O + 4 Fe(OH)3 + 8 H+ -> CO2 + 4 Fe2+ + 11 H2O, using the `g2` OM pool and amorphous Fe(OH)3 (`amfe3`).
  Rate = vmax * g2/(g2+Kg2) * amfe3/(amfe3+Kfe) * O2 inhibition * NO3- inhibition * min(2 sat, 1),
  with an analytic Jacobian. vmax and Kfe are `krxn1_ext_all` / `krxn2_ext_all` for `g2fe3` and are
  provisional (200 mol CH2O m-3 yr-1 and 100 mol m-3); they need calibrating against paddy incubations.
  Proton consumption is not written explicitly: SCEPTER gets pH from charge balance.
  Enable it by listing `g2fe3` in `extrxns.in`.
- Fe2+ on the cation-exchange complex. Fe2+ is treated like Mg2+ in the general-charge exchange solver
  (`get_maqads_all_v4`/`v4a`), with log KH\Fe = log KH\Mg + 0.158, from the Appelo and Postma (2005)
  selectivities (phreeqc.dat FeX2 0.44 vs MgX2 0.6). `cec.in` keeps its five columns; the Fe2+ value
  follows whatever Mg value is read there. Exchange fluxes are written to `flx_aqex-fe2.txt`.
- Siderite (`sid`, FeCO3) as a new solid: log K -10.89 and delta H -10.38 kJ/mol (phreeqc.dat),
  molar volume 29.378 cm3/mol, neutral-mechanism kinetics from Palandri and Kharaka (2004)
  (log k -8.90 mol m-2 s-1, Ea 62.76 kJ/mol). It is allowed to precipitate in `data/2ndslds_paddy.in`
  and must also be listed as a tracked solid.

- Sulfate reduction `g2so4` (2 CH2O + SO42- + H+ -> 2 CO2 + HS-): Monod in OM and SO42-, O2 and NO3-
  inhibition as `g2fe3`. vmax 100 mol CH2O m-3 yr-1 and K(SO42-) 1e-4 mol L-1 (`krxn1`/`krxn2`), provisional.
- New gas `pch4`: Henry constant 1.4e-3 mol L-1 atm-1 at 25 C with 1600 K temperature dependence (Sander 2015),
  aqueous diffusion 1.49e-9 m2 s-1 at 15 C (Jahne et al. 1987), gas-phase diffusion as CO2.
  The Python drivers now write pch4 = 1.9e-6 atm to `atm.in` by default.
- Methanogenesis `g2ch4` (2 CH2O -> CH4 + CO2, acetoclastic and hydrogenotrophic lumped): Monod in OM,
  inhibited by O2, NO3-, SO42- (constant = K(SO42-) of `g2so4`) and amorphous Fe(III) (`krxn2`, 30 mol m-3),
  so it follows the redox sequence. vmax 100 mol CH2O m-3 yr-1, provisional.
- Aerobic CH4 oxidation `ch4o2` (CH4 + 2 O2 -> CO2 + 2 H2O): dual Monod in CH4(aq) (K 5e-6 mol L-1)
  and O2(aq) (K 2e-5 mol L-1), vmax 500 mol m-3 yr-1, provisional.
- Hydroxy green rust (`gr`, Fe(II)Fe(III)2(OH)8) as a new solid, formed by the extra reaction `fe2gr`:
  Fe2+ + 2 Fe(OH)3 (amfe3) -> Fe3(OH)8 + 2 H+, log K -9.98 (from Bourrie et al. 1999, log K(gr) 45.8 and
  log K(Fe(OH)3) 13.02 per Fe). Rate = k * ln(omega) * f, with omega = K [Fe2+]/[H+]2, f = amfe3/(amfe3+K)
  when forming and gr/(gr+K) when dissolving (k 1e4 mol m-3 yr-1, K 1 mol m-3, provisional). Being linear in
  ln(omega) it acts as the Fe(OH)3-Fe2+ pH-stat seen in flooded soils (Ponnamperuma 1972): pH tends to
  4.99 - 0.5 log[Fe2+]. Its mineral dissolution/precipitation rate constant is set to zero, because as an
  ordinary mineral it would need aqueous Fe3+ (~1e-20 M at pH 7) and either stalled the solver or never formed.
  List `gr` as a tracked solid and `fe2gr` in `extrxns.in`.
- Plant-mediated CH4 transport `ch4pl`: first order in dissolved CH4 above atmospheric equilibrium,
  R = k(z) * poro * sat * 1e3 * KH * (pCH4 - 1.9e-6), with k(z) = k0 exp(-z/zroot) following root density
  (k0 50 yr-1, zroot 0.1 m, provisional; the depth profile is applied in the main program after each
  `coefs_v2` call). Half of the CH4 taken up is oxidised to CO2 in the rhizosphere (Pox = 0.5,
  Walter and Heimann 2000) and half is emitted, so CH4 emitted via plants = 0.5 x the `ch4pl` flux.
- CH4 ebullition `ch4eb`: R = k * poro * sat * 1e3 * KH * max(pCH4 - pthr, 0), k 100 yr-1, pthr 0.4 atm
  (about the 500 uM bubble threshold of Walter and Heimann 2000), provisional. Ebullited CH4 leaves directly.
- Iron monosulfide (`fes`, molar volume 20.4 cm3/mol) formed by `fe2fs`: Fe2+ + HS- -> FeS + H+,
  log K 3.915 (FeS(ppt), phreeqc.dat). Rate = k * ln(omega) with omega = K [Fe2+][HS-]/[H+], limited by FeS
  present when dissolving (k 1e3 mol m-3 yr-1, K 1 mol m-3, provisional). FeS is oxidised by `fesox`:
  FeS + 2 O2 -> Fe2+ + SO42-, rate = k [FeS] [O2(aq)] with k 1e6 L mol-1 yr-1 (about a day at air
  saturation, provisional). Like `gr`, FeS has no mineral rate law of its own. Track `fes` as a solid
  (added automatically with `--ch4`) and list `fe2fs,fesox` in `extrxns.in`.
- N2O reduction for paddies `g2n2p` (2 N2O + CH2O -> 2 N2 + CO2): Monod in OM and N2O(aq) (K 5 uM, the
  micromolar affinity of soil N2O reductase), inhibited by O2 (2.52e-5 M), by NO3- (1 mM: the N2O/(N2O+N2)
  ratio rises mainly at millimolar nitrate) and by low pH through 1/(1 + (H+/H50)^2) with pH50 = 6.0 (`krxn2`).
  vmax `krxn1` (default as Maggi et al. 2008; calibrated 300 mol m-3 yr-1, below). Listed with `amo2o,g2n21`.
- Plant-mediated N2O transport `n2opl`: as `ch4pl` (same root-density profile), with the CH4 conductance scaled
  by the ratio of gas diffusivities in air (`krxn1` = 0.76) and no rhizosphere loss.
- Urea (`urea`, CO(NH2)2 + H2O -> 2 NH3 + CO2, NH3/NH4+ carried as `pnh3`) as a new solid that decays with a
  5-day turnover at 25 C (Q10 2). Supplied with the OM input; driver option `--nfert=<kg N/ha>`.
- Floodwater gas exchange. SCEPTER closes the top boundary to dissolved gases and lets gases leave only
  through soil air, so a fully saturated column had no surface exchange at all (no O2 in, no CO2, CH4 or
  N2O out except by leaching). `gasa_diff` now adds a transfer velocity across standing water, read from
  `floodwater.in` (header line, then the value in m/yr; default 0 keeps the upstream behaviour). It is
  weighted by (sat(1) - 0.9)/0.1, so it acts only while the surface is flooded. The driver writes 36.5 m/yr
  (0.1 m/d, about the diffusive boundary layer at the soil surface; `--kfw`).
- Run-time rate constants: `krxn_ext.in` in the run directory (header line, then `name krxn1 krxn2`; -999 keeps
  the default) overrides `krxn1_ext_all`/`krxn2_ext_all` at the end of `coefs_v2`. Driver option `--k=<file>`.
  `g2n21`/`g2n0` now take their vmax from `krxn1_ext_all` too (same default).
- Bug fix in `calc_rxn_ext_dev_3`: `maqft_loc` for NO3- is the total concentration, not a total/free ratio
  (as for SO42-), so `no3x = maqx*maqft` squared nitrate and made denitrification and every NO3- inhibition
  term ~1e4 times too weak. `no3x` is now `maqx_loc`.
- `makefile_paddy`: as `makefile_ex11` but without `-Ddisp_cnst`. That flag (needed for the PHREEQC ex11
  benchmark) gives every aqueous and gas species the same diffusivity, 0.066 m2/yr, which made gas-phase
  diffusion ~1e4 times too slow: drained soil stayed anoxic and made methane. Paddy runs use `scepter_paddy`.

## Tests

The column tests below up to the AWD sensitivity used the ex11 build (constant diffusivity for all species,
see `makefile_paddy`), so their numbers are superseded by the calibration runs at the end.

- `compare_ex11_phreeqc.py`: checks the cation-exchange benchmark against PHREEQC example 11
  (via phreeqpython). Outlet Na, K, Ca, Cl agree with PHREEQC to within about 0.06 mM (RMSE under 0.015 mM).
  Repeated runs of the unmodified upstream code differ from each other by up to about 0.015 mM,
  so bit-for-bit comparison between builds is not meaningful.
- `run_flooded_column.py`: 0.3 m water-saturated column with 3 wt% amorphous Fe(OH)3 and OM input.
  With `g2fe3` on, about 17 mol m-2 of Fe(III) is reduced in 0.5 yr below the oxic surface layer,
  and the mass balance closes (residual ~1e-11).
  With `--fe2buf` (CEC 10 cmol/kg on `inrt`, Na/K/Mg added, siderite tracked), pH stays at 6.5-7.8
  instead of ~12: Fe2+ takes 25-37 % of the exchange sites and displaces Ca and Mg into solution.
  Fe2+ mass balance closes (residual ~2e-9 mol m-2 yr-1).
  With `--ch4` (adds SO42-/HS-, CH4 and 0.2 mM CaSO4 irrigation water) and
  `fe2o2,g2fe3,g2so4,g2ch4,ch4o2`: sulfate falls from 0.2 mM to <1e-6 mM by 15 cm, CH4 builds to ~0.1 atm
  equivalent at depth, ~0.15 mol CH4 m-2 yr-1 is produced and ~12 % of it is oxidised in the top 5 cm.
  Mass balances close for CH4, SO42-, Fe2+ and OM.
  Adding `fe2gr` (`fe2o2,g2fe3,g2so4,g2ch4,ch4o2,fe2gr`) keeps pH at 6.56-6.89 for the whole 0.5 yr run
  (without it pH reaches ~11.7), with 3-22 mol m-3 green rust below the surface layer. CH4 production
  (0.15 mol m-2 yr-1, ~17 % oxidised) is barely changed. Run time ~4 min.
  Adding `ch4pl,ch4eb`: of 0.152 mol m-2 yr-1 produced, 0.065 goes into plant transport (half emitted),
  0.045 leaves through the bottom of the column with the 0.5 m yr-1 infiltration (the `adv` column of
  `flx_gas-pch4.txt`; an artefact of the shallow test column), 0.023 is oxidised aerobically and 0.020 is still
  accumulating. Diffusion to the air through the saturated surface is negligible, so in this column all
  emitted CH4 goes through plants (field: ~60-90 %, the rest by ebullition and diffusion). pCH4 peaks at 0.06 atm, so ebullition does not start. With the
  threshold lowered to 0.03 atm as a code test, pCH4 is capped at 0.034 atm and ebullition takes 0.048.
  CH4 balances close (residual ~1e-17); pH stays at 6.56-6.84.
  Adding `fe2fs,fesox`: all 0.14 mol m-2 yr-1 of sulfide from sulfate reduction precipitates as FeS
  (0.06-0.46 mol m-3), free HS- falls from ~5e-5 M to ~5e-8 M, and 0.037 mol m-2 yr-1 of FeS is re-oxidised
  near the surface, recycling sulfate. S balances close (residual ~4e-14); pH 6.55-6.79; run time ~6.5 min.
- Flood/drain cycling (`--awd`): uses SCEPTER's seasonal forcing (`season` on, `Wet_temp.in` = surface
  saturation, water table at 0.5 m). Flooded for 0.08 yr, then ~1 week drains to surface saturation 0.5 and
  ~2 week refloods until 0.5 yr. `--season` runs the same machinery kept flooded, as the control.
  Season-integrated (0-0.5 yr, all reactions on), control vs AWD, mol m-2:
  CH4 produced 0.035 vs 0.029 (-16 %); emitted (half of plant transport; diffusion and ebullition ~0)
  0.0069 vs 0.0042 (-39 %, in the 30-50 % range reported for AWD); aerobic CH4 oxidation 0.0049 vs 0.0057;
  Fe(III) reduced 7.6 vs 6.9; FeS re-oxidised 0.018 vs 0.046, which recycles sulfate so sulfate reduction rises
  (0.079 vs 0.100). pH 6.55-6.87 vs 6.56-7.00. Budgets close. Run time ~14 min each. To integrate, sum
  each row of the `flx_*.txt` files times the time since the previous row.
- CH4 magnitude sensitivity (`--fe3=<wt fraction>`, `--om=<g C m-2 yr-1>`, flooded 0.5 yr, all reactions):
  season CH4 emission 0.11 g m-2 (base: 3 wt% amfe3, 300 g C), 0.36 (1 wt%), 0.28 (600 g C), 1.0 (1 wt% + 600 g C),
  4.5 g m-2 (0.3 wt% + 600 g C; ebullition starts). Field: ~10-40 g m-2 per season. Even in the last case only
  ~10 % of the OM added is respired, so the remaining gap points to the provisional methanogenesis/OM rate
  constants (vmax, K_OM) rather than a missing process.

- Calibration against field data (`calib_metrics.py <run dir> [days] [kg N/ha]`). Protocol: a 3-yr drained
  fallow spin-up (`--fallow`, surface saturation 0.5, water table 0.5 m; aerobic reactions only, plus `amo2o` for
  N runs) builds the labile OM stock (g2 80-175 mol m-3) and leaves the soil oxic; then a 113-day season
  (IPCC 2019 default length) restarts from it (`--rstrt=<tag>`), with 0.5 wt% amorphous Fe(OH)3 and
  300 g C m-2 yr-1 OM input. Fitted constants (`--k=paddy/krxn_calibrated.in`): vmax 300 mol CH2O m-3 yr-1 for g2fe3, g2so4, g2ch4
  (and g2n21, g2n2p in N runs), Fe(III) inhibition of methanogenesis 5 mol m-3, plant transport k0 300 yr-1.
  Safe AWD (`--awd`): drains to a water table at 15 cm and surface saturation 0.6 (`--zdry`, `--sdry`).

  | Metric | Flooded | AWD | Field target |
  |---|---|---|---|
  | CH4 emission, kg CH4 ha-1 d-1 | 1.17 | 0.64 | 1.19 (IPCC 2019 default) |
  | AWD / flooded CH4 | | 0.55 | 0.55 (IPCC multiple drainage) |
  | Plant share of CH4 emission | 70 % | 63 % | 60-90 % |
  | CH4 onset (10 % of peak rate) | 21 d | 16 d | 1-3 weeks |
  | Peak porewater CH4 / Fe2+, mM | 0.82 / 1.2 | 0.86 / 1.1 | 0.1-1 / 0.5-3 |
  | N2O EF (120 kg N/ha urea, minus no-fertiliser run) | 0.37 % | 2.2 % | 0.3 % / 0.5 % (IPCC) |
  | N2O background, kg N/ha/season | 0.95 | 3.3 | ~0.2-1 |

  With N on, CH4 is 1.14 kg/ha/d flooded and 0.58 AWD (ratio 0.51). 84-97 % of N2O made by denitrification is
  reduced to N2 when flooded. About 31 % of urea N volatilises as NH3 (field 10-40 %).

- Recalibration on the reference paddy soil (2026-10-24; `krxn_calibrated.in` now holds this set). Reference soil
  (`--refsoil`, CEC 15, 2 wt% goethite `--gt=0.02`, green rust with 2 % bioavailable Fe(III), rice root respiration
  30 mol C m-2 yr-1, no siderite nucleation), spin-ups `sgt12_450` (1.2 wt% amfe3, 450 g C) and `sgt08_300`
  (0.8 wt%, 300 g C), 120 kg N/ha urea. Only the Fe(III) inhibition constant of methanogenesis (g2ch4 krxn2) was
  fitted; season CH4 rises almost linearly with it:

  | g2ch4 Ki (mol m-3) | 5 | 6.5 | 7 | 8 | 10 | 13 | 16 |
  |---|---|---|---|---|---|---|---|
  | CH4, kg/ha/d, 1.2 % Fe soil | 1.04 | 1.18* | 1.24* | 1.44 | 1.66 | 1.93 | 2.16 |
  | CH4, kg/ha/d, 0.8 % Fe soil | | 1.11* | 1.16* | 1.34 | 1.53 | 1.77 | |

  (* 15-cell grid, `--nz=15`, with its own 15-cell spin-ups `sgt12_450n15`, `sgt08_300n15`; this grid gives ~6 %
  less CH4 than 30 cells at Ki 10 and runs ~3x faster, ~50 min per season, so it is the production grid.)
  With aerobic respiration of labile OM (`g2o2`) listed, Ki 6.5 gives 1.24 and 1.17 kg/ha/d (mean 1.21; Ki 6: 1.18, 1.11),
  so `krxn_calibrated.in` uses Ki 6.5 and production runs list `g2o2`. Without g2o2: at Ki 7 the two soils average 1.20 kg/ha/d (IPCC 1.19), with pH 6.5-6.6, porewater Fe2+ 0.6-2.5 mM and pCO2
  0.07-0.27 atm through the season, plant share 61-67 %. Still off: CH4 onset 34-50 d (field 1-3 weeks) and
  the flooded N2O EF 0.4-0.7 % (IPCC 0.3 %).

## Known gaps (next steps)

- The green rust pH-stat makes pH depend on green rust thermodynamics; its log K is uncertain by about
  a unit, which moves the buffered pH by ~0.5. Should be checked against paddy porewater pH-Fe2+ data.
- The calibration fits one site-type (0.5 wt% Fe(III), 300 g C) to global-default IPCC factors; it needs
  testing against site time series (CH4 flux, porewater Fe2+/CH4) and the ERW-in-rice datasets.
- N2O under AWD is ~4x too high (EF 2.2 % vs 0.5 %) and the flooded background is at the top of the field
  range. Likely causes: no nitrifier N2O pathway or NO/N2O partitioning, N2O reductase recovery after reflooding,
  and the nitrate built up in the fallow spin-up.
- The fallow spin-up is not at steady state (3 yr) and has no crop; with the N species it accumulates 1-1.5 mM
  nitrate and acidifies to pH ~5.9.
- No rice growth curve or root O2 release for plant transport; k0 and zroot need setting against field data.

- Siderite does not nucleate: SCEPTER's default precipitation rate scales with the amount of mineral
  already present (seed 1e-20 mol m-3), so it stays at the seed even at log omega ~2.7. Paddy
  porewaters are often siderite-supersaturated, so this may be acceptable, but a nucleation or
  seeded-precipitation option is worth testing against data.
- Without exchange (no `--fe2buf`) the column still rises to pH ~12, so CEC must be on for paddy runs.
- No pyrite formation from FeS (slow, years); FeS is the only sulfide sink.
- Draining is represented through the saturation profile and the floodwater transfer velocity: there is no
  explicit floodwater layer or water-table dynamics, and no change in infiltration between flooded and drained periods.
