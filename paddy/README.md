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

## Known gaps (next steps)

- Siderite does not nucleate: SCEPTER's default precipitation rate scales with the amount of mineral
  already present (seed 1e-20 mol m-3), so it stays at the seed even at log omega ~2.7. Paddy
  porewaters are often siderite-supersaturated, so this may be acceptable, but a nucleation or
  seeded-precipitation option is worth testing against data. Green rust and FeS are not defined.
- Without exchange (no `--fe2buf`) the column still rises to pH ~12, so CEC must be on for paddy runs.
- No sulfate reduction, methanogenesis, CH4 gas or CH4 oxidation yet.
- Flood and drain cycling is not yet driven through the seasonal wetness forcing (`Wet_temp.in`).
