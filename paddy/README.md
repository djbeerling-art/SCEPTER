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

## Tests

- `compare_ex11_phreeqc.py`: checks the cation-exchange benchmark against PHREEQC example 11
  (via phreeqpython). Outlet Na, K, Ca, Cl agree with PHREEQC to within about 0.05 mM.
  Repeated runs of the unmodified upstream code differ from each other by up to about 0.015 mM,
  so bit-for-bit comparison between builds is not meaningful.
- `run_flooded_column.py`: 0.3 m water-saturated column with 3 wt% amorphous Fe(OH)3 and OM input.
  With `g2fe3` on, about 17 mol m-2 of Fe(III) is reduced in 0.5 yr below the oxic surface layer,
  and the mass balance closes (residual ~1e-11).

## Known gaps (next steps)

- Porewater pH rises to ~12 in the test column because nothing yet buffers the alkalinity from
  Fe(III) reduction: Fe2+ is not on the exchange complex, siderite did not precipitate within the
  run, and Fe(OH)2 / green rust are not defined. Real paddies converge to pH ~6.5-7.
- No sulfate reduction, methanogenesis, CH4 gas or CH4 oxidation yet.
- Flood and drain cycling is not yet driven through the seasonal wetness forcing (`Wet_temp.in`).
