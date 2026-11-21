# Verification (GMD Fig. 3)

Run on unpruned one-year runs (season + fallow) at archetype irr061, flooded and AWD, 0 and 20 t/ha
(`paddy/multiyear/run_chain.sh` without its pruning steps, restarting from a settled state).

- `ex11_data.py`: PHREEQC example 11 cation-exchange benchmark (as `paddy/compare_ex11_phreeqc.py`), writes `ex11.json`.
  Usage: `SCEPTER_OUT=<ex11 outdir> python3 paddy/verification/ex11_data.py` (edit the json path at the end).
- `phreeqc_respeciate.py <out.json> <prof dir>:<index>:<tag> ...`: re-speciates every depth of the model porewater
  (`charge_balance-<index>.txt`) in PHREEQC (phreeqc.dat, pH by charge balance) and compares pH and saturation of
  calcite, siderite, Fe(OH)3(a) and SiO2(a) with `sat_sld-<index>.txt`. `TEMP_C` sets the temperature (run value).
- `mass_balance.py <out.json> <run dir> ...`: C, Fe, N, Ca, Mg budgets from `int_flx_*`; checks that internal
  reaction terms cancel and that the storage term matches the change in profile inventories.
- `fig3_verification.py`: draws the figure from the three json files in its own directory.

Results (v16, 28 Oct 2026): ex11 RMSE <= 0.019 mM; pH within 0.04 (105 porewaters); mean SI offsets (ideal
solution) calcite +0.13, Fe(OH)3(a) -0.06, SiO2(a) 0.00, siderite +0.08 after the model's +1 log-unit apparent
solubility; Fe, N, Ca, Mg reaction closure <= 1e-11 of gross flux, storage <= 1e-7. C leaks 4% of aerobic OM
decomposition (upstream oxalate/acetate released by `g2` but not tracked when those species are off).
