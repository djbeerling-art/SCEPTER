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

## Tests

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
  0.045 diffuses out of the top, 0.023 is oxidised aerobically and 0.020 is still accumulating; plants carry
  ~42 % of the 0.077 mol m-2 yr-1 emitted. pCH4 peaks at 0.06 atm, so ebullition does not start. With the
  threshold lowered to 0.03 atm as a code test, pCH4 is capped at 0.034 atm and ebullition takes 0.048.
  CH4 balances close (residual ~1e-17); pH stays at 6.56-6.84.

## Known gaps (next steps)

- The green rust pH-stat makes pH depend on green rust thermodynamics; its log K is uncertain by about
  a unit, which moves the buffered pH by ~0.5. Should be checked against paddy porewater pH-Fe2+ data.
- Plant transport carries only ~42 % of emitted CH4 in the test column, against the 60-90 % usually measured.
  The column has no standing floodwater (which slows diffusive loss) and no rice growth curve, so k0, zroot and
  a seasonal plant-growth factor need setting against field data. Root O2 release is not yet represented.

- Siderite does not nucleate: SCEPTER's default precipitation rate scales with the amount of mineral
  already present (seed 1e-20 mol m-3), so it stays at the seed even at log omega ~2.7. Paddy
  porewaters are often siderite-supersaturated, so this may be acceptable, but a nucleation or
  seeded-precipitation option is worth testing against data. FeS is not defined.
- Without exchange (no `--fe2buf`) the column still rises to pH ~12, so CEC must be on for paddy runs.
- No FeS / pyrite formation from HS-, so HS- accumulates.
- Flood and drain cycling is not yet driven through the seasonal wetness forcing (`Wet_temp.in`).
