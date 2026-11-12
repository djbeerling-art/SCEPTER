# Arsenic in Riceland SCEPTER: Step 1 (arsenic branch)

This adds trace arsenic to the paddy model. It is kept on a separate branch so the v1.0 GMD code is untouched. Runs
that do not list the arsenic species or reactions give bit-identical output to v1.0. That was checked on a flooded
column, apart from the two extra charge-balance columns in `ph.txt`.

## What is tracked

| Name | Kind | Meaning |
|---|---|---|
| `as3`, `as5` | solutes | dissolved As(III) (H3AsO3) and As(V). Both are trace: no charge, no aqueous complexes. Diffusivity is that of H4SiO4 for As(III) and of H2AsO4- for As(V) (Li and Gregory 1974). |
| `as3ad`, `as5ad` | solids | As(III) and As(V) sorbed on Fe(III) oxide. Molar weight 74.92 (As), so wt fractions are As contents. They react only through the extra reactions below. |
| `as3sp`, `as5sp` | reactions | sorption, relaxing towards the Si-competitive equilibrium (below) |
| `as5rd`, `as5ra` | reactions | microbial As(V) reduction, of dissolved and sorbed As(V). First order, O2-inhibited as `g2fe3`. |
| `as3ox`, `as3oa` | reactions | As(III) oxidation by O2, of dissolved and sorbed As(III). First order, Monod in O2(aq). |
| `as3up` | reaction | rice uptake of As(III) through Lsi1, competitive with Si (Ki = Km of `siup`). Root profile as `siup`. Placeholder constants. |

## Sorption equation (`as3sp`, `as5sp`)

```
sorbed_eq [mol m-3] = c [mol L-1] * K * ( f_site * 0.2 * amfe3 / (1 + KSi * cSi^n)  +  Sx )
log K   = a + b (pH-7) + c (pH-7)^2        log KSi = aS + bS (pH-7)
rate    = krxn1 * min(2 sat, 1) * (sorbed_eq - sorbed)
```

- 0.2 is the number of weak sites per mol Fe (Dzombak and Morel 1990).
- `f_site` comes from `as_sites.in` and defaults to 1. It scales the site density for aged or P- or OM-coated oxide.
- `Sx` (krxn2) is the Si-insensitive sites in mol m-3, standing for crystalline Fe, clays and OM. It defaults to 0 and is soil-specific.
- Sorbed As follows `amfe3`, so it is released as the ferrihydrite is reductively dissolved (`g2fe3`).
- The constants were fitted to the PHREEQC diffuse-layer model in `paddy/arsenic/` (`grid_phreeqc.py`, `fit_kd.py`). That model uses USGS `wateq4f.dat` HFO constants for arsenite, arsenate and phosphate, and Swedlund and Webster for Si. The grid has 882 points: pH 5–8, Si 0.02–2 mM, As 2–200 µM, ferrihydrite Fe 0.02–0.23 mol per L porewater, and P 0.01 mM. Fitted constants and fit quality:

| | a | b | c | aS | bS | n | RMSE (log10) | max error (log10) |
|---|---|---|---|---|---|---|---|---|
| As(III) | 5.1882 | -0.0224 | -0.0171 | 3.8028 | 0.0104 | 0.8685 | 0.012 | 0.058 |
| As(V) | 7.2189 | -1.0735 | -0.2022 | 9.3859 | 0.1815 | 1.9115 | 0.052 | 0.223 |

pH and Si are clamped to the fitted range. In the fitted model:

- Raising Si from 0.3 to 1.5 mM raises dissolved As(III) 3.57 times. The full DLM gives 3.6.
- Losing 80% of the ferrihydrite raises it 5 times.
- As(V) rises 54 times from pH 6 to 7.5.

## Assumptions to revisit

- As is trace relative to sites, so As(III) and As(V) do not compete with each other. Loading in the fit grid is below 5%.
- Si sorption is not removed from solution. Porewater Si comes from SCEPTER's own balance.
- SCEPTER carries no dynamic P. Phosphate enters only as a fixed background competitor, through As sorption constants refitted at labile P/Fe 0.02 and 0.1 (`kd_constants_Pr*.json`, `--askd`). There is no carbonate or DOC competition.
- Scope: the model predicts changes in root As(III) uptake under flooded conditions. It does not yet predict total root arsenic uptake or grain inorganic arsenic. Before extending to alternate wetting and drying, add a simple As(V) uptake pathway with phosphate inhibition (As(V) matters more during drainage).
- Not represented (state as limitations): Mn oxides (As(III) oxidation and sorption in drained phases), Ca/Mg-arsenate precipitation and Ca-bridged arsenate sorption, As(V) root uptake via phosphate transporters, transient P release with Fe reduction, and grain partitioning.
- The redox rate constants are provisional: `as5rd` 50 /yr and `as3ox` 1000 /yr with O2 half-saturation 2e-5 M. They are not calibrated.
- The `as3up` constants are placeholders: Vmax 5e-3 mol m-2 yr-1 and Km 30 µM.

## Running

```bash
python3 paddy/run_flooded_column.py <tag> <reactions> <yr> ... --as=<mg/kg>[:<ug/L irrigation>[:<As(V) fraction>]] \
    [--asx=<mol/m3>] [--assite=<fraction>]
```

- `--as` adds the species and the sorption and redox reactions. It puts the labile soil As in as sorbed As(V), the state of oxidised soil before flooding.
- To get uptake, list `as3up` together with `siup`.
- `paddy/arsenic/as_profile.py <run dir> <z>` prints pH, Si, As(III), As(V), amfe3 and sorbed As through time at depth z.

## First test: one flooded season, 0.3 m column, no spin-up (uncalibrated)

Setup: reference soil with 3 wt% amfe3 and 6 mg/kg labile As. Reactions: `fe2o2,g2fe3,g2gr,g2so4,g2ch4,ch4o2,fe2gr,ch4eb,fe2fs,fesox,g2o2,groxy,sidox`. Options: `--ch4 --n --kaol=0.05 --refsoil --gt=0.02 --season --k=paddy/krxn_calibrated_v25_cal.in --nz=15 --ztot=0.3`. Results at 10 cm after 0.3 yr:

| Run | pH | Si (µM) | As(III) (µM) | As(V) peak (µM) |
|---|---|---|---|---|
| control | 7.27 | 266 | 0.053 | 0.11 (day 7) |
| basalt 20 t/ha | 8.31 | 1333 | 0.210 (4.0×) | 2.9 (weeks 2–6) |
| 1 wt% amfe3 | 7.28 | 265 | 0.177 (3.3×) | 0.32 |

- Sorbed As(V) is reduced to As(III) over about six weeks after flooding.
- Si competition alone accounts for 3.5× of the 4.0× As(III) rise with basalt.
- The basalt pH rise also mobilises As(V) for a few weeks, before it is reduced.
- Column mass balance closes. Residuals are below 1e-11 mol m-2.
- Major-element chemistry with and without As differs by about 1%. That is solver path; As is trace.
- Plugging the basalt run into the competitive uptake ratio (Si Km 300 µM) gives 4.0 × (1 + 266/300)/(1 + 1333/300) = 1.39. That is a 39% rise in As(III) uptake.
- This is a smoke test, not a prediction. There was no spin-up, the basalt pH is above the fitted range, and the redox constants are uncalibrated.

## Calibration and response-map tools (2 Nov 2026)

- **Si-insensitive sites on goethite.** `as3sp` and `as5sp` add f_gt sites per mol goethite Fe (`gt`) to krxn2. f_gt defaults to 0.02 and is the second value in `as_sites.in`, set with `--assite=<f>:<fgt>`. Crystalline Fe (dithionite minus oxalate) is therefore the soil input for the Si-insensitive pool.
- **Pot tests.** Seyfferth & Fendorf 2012 Richvale and Hawi pots, run at the SI values (oxalate Fe, DCB Fe, As spike), give absolute porewater As(III) within about 3× without fitting. The digitised Fig. 1 is in `paddy/arsenic/seyfferth2012_fig1_digitised.csv`. Two problems remain: silica gel dissolves too slowly unless amsi is scaled (production uses `--kscale=...,amsi:30`), and `--refsoil` irrigation gives too high a pH for these pots.
- **Threshold test** (`paddy/arsenic/threshold_table.py`). It compares the porewater As ratio with the Si-inhibition threshold T = (1+Si_trt/Ki)/(1+Si_ctrl/Ki) at Ki = 300 µM. The grain As direction is right in 12 of 15 treatments. Steel slag (He 2020), the basalt analogue, lies above the threshold. Richvale silica gel lies on the threshold but grain As fell, which points to Si effects beyond competition.
- **Response-map sweep** (`paddy/arsenic/run_as_sweep.sh`). 27 soils (ferrihydrite × goethite × irrigation alkalinity) each get a spin-up, a settle, then control and 20 t/ha basalt chains. `as3up` runs in the season only (`XRSS` in `run_chain.sh`). `paddy/arsenic/as_sweep_summary.py` reports pH, Si, As(III), R_As, T, the competitive and simulated uptake ratios, and the class (suppress, mobilise or uncertain). Grain iAs is left empty until translocation is modelled.
- **Smoke test.** One soil (1.2% ferrihydrite, 2% goethite, reference irrigation; 1 settle year, 1 chain year) ran end to end. Basalt raised Si 7× but As(III) only 1.27×, against a threshold of 4.2. The uptake ratio was 0.26 simulated and 0.30 competitive, so the class is suppress. The goethite pool dominates As retention in this soil.

## Goethite x silica test (2 Nov 2026)

Fresh 0.3 m columns, 1.2 wt% ferrihydrite, 10 mg/kg As, 100 ug/L irrigation As, 0.12 yr flooded. Goethite 0, 2 and
6 wt%, each without and with silica gel (`--amsi=0.077 --kscale=amsi:30`), so Si rises about 7x (220 to 1570 uM) with
no pH change. All six runs start from the same porewater As (0.80 uM As(III), 0.09 uM As(V)). Root zone 0-0.2 m,
time means over the run (`paddy/arsenic/gt_si_partition.py`).

| Goethite | As(III) control (uM) | As(III) +Si (uM) | Absolute change (uM) | R_As | R_As predicted from Kd | Ferrihydrite share of As(III) capacity |
|---|---|---|---|---|---|---|
| 0 | 0.081 | 0.388 | +0.307 | 4.78 | 4.86 | 0.995 |
| 2% | 0.045 | 0.074 | +0.029 | 1.64 | 1.70 | 0.51 |
| 6% | 0.025 | 0.031 | +0.006 | 1.22 | 1.27 | 0.26 |

- Si cuts the ferrihydrite Kd for As(III) about 5x (1230 to 250). The goethite Kd changes by less than 1%. Goethite is a genuinely Si-insensitive pool.
- R_As is set by the ferrihydrite share of sorption capacity, phi: R_As = 1/(1 - phi(1 - 1/D)), D about 5. Total Fe is not the control.
- Sorbed As(III) moves from ferrihydrite to goethite sites with Si (end of run, mmol m-2): at 2% goethite 14.9/14.3 to 5.0/24.0; at 6% 7.6/22.0 to 1.9/27.5; without goethite 28.8/0.1 to 28.1/0.7.
- That split is diagnostic. The model holds one sorbed pool per As species at the competitive equilibrium of both site types, so the transfer is instantaneous by construction. A kinetic or irreversible transfer needs separate tracked pools.
- Root-zone pH is about 7.6 in these runs (fresh `--refsoil` columns), inside the fitted range. A first version of the script averaged the 0-0.2 m mask over the whole column and reported pH 5.05 and Si and As(III) two-thirds too low; ratios and inventories were unaffected.

## Ferrihydrite-capacity predictor and silica test from settled soils (2 Nov 2026)

Predictor (David Beerling): if Si cuts the ferrihydrite Kd by q and goethite is unaffected, R_As = 1/((1-f) + f/q),
f = ferrihydrite share of sorption capacity in the control. It assumes linear equilibrium sorption and a conserved
exchangeable As inventory. `as_sweep_summary.py` now reports f_fh, q, q_Si, R_As_fq, R_As_cap and R_inv (exchangeable
As(III) inventory, dissolved + sorbed, B/C) from the model's own time-resolved Kd values.

Test: `paddy/arsenic/run_si_test.sh` (spin-up, 2 settle years, then one season + fallow from the settled state with and
without 7.7 wt% silica gel in the plough layer; `add_silica.py`). 1.2% ferrihydrite; goethite 0.5, 2 and 6%;
irrigation alkalinity 1 and 5, and alkalinity 5 with 2% soil calcite (`SOILX=--soilcc=0.02`). Silica raised Si about
4.9x (240 to 1170 uM) with no pH change. Summary: `python3 paddy/arsenic/as_sweep_summary.py <out>/paddy si`.

| Soil | Season pH | f | q | R_As (model) | R_As f/q rule | R_inv | log10(R_As / rule) |
|---|---|---|---|---|---|---|---|
| gt 0.5%, alk 1 | 6.44 | 0.86 | 3.87 | 1.92 | 2.73 | 0.83 | -0.15 |
| gt 2%, alk 1 | 6.46 | 0.60 | 3.84 | 1.29 | 1.80 | 0.84 | -0.15 |
| gt 6%, alk 1 | 6.46 | 0.32 | 3.81 | 1.01 | 1.31 | 0.84 | -0.11 |
| gt 0.5%, alk 5 | 6.57 | 0.86 | 3.84 | 1.88 | 2.74 | 0.84 | -0.16 |
| gt 2%, alk 5 | 6.59 | 0.60 | 3.82 | 1.27 | 1.80 | 0.85 | -0.15 |
| gt 6%, alk 5 | 6.59 | 0.32 | 3.79 | 1.01 | 1.31 | 0.85 | -0.11 |
| gt 0.5%, alk 5, 2% calcite | 6.70 | 0.86 | 3.78 | 1.85 | 2.74 | 0.84 | -0.17 |
| gt 2%, alk 5, 2% calcite | 6.71 | 0.60 | 3.75 | 1.25 | 1.78 | 0.85 | -0.15 |

(The calcareous 6% goethite soil failed in its second settle year and was not rerun.)

- The rule orders the soils correctly and carries the goethite effect, but over-predicts the seasonal R_As by 30-45%
  (0.11-0.17 log units, least where f is small).
- About half of the gap is the inventory assumption: with Si the exchangeable As(III) inventory in the root zone falls
  15-17% in every soil, because desorbed dissolved As(III) is oxidised (`as3ox` +78% in the 2% goethite soil). The rest is
  depth and time variation in Kd and Si (silica dissolves near the surface): per date, R_inv x rule still exceeds R_As
  by 7-20%.
- pH: settled flooded seasons converge to pH 6.4-6.7 whatever the irrigation alkalinity or 2% calcite, and the results
  barely move over that range (q 3.87 to 3.75, R_As within 4%). The fresh columns (pH 7.6) gave the same dependence on f.
  Settled seasons at pH 7-7.5 were not reached with these settings.
- As(V): dissolved As(V) is about 8% of As(III) over the season and Si raises it about 10%, so the pH sensitivity of
  As(V) sorption has little leverage on exposure here.
- Silica alone raises Si about 5x with no alkalinisation; T = 2.74 exceeds R_As in every soil, so all are 'suppress'
  (competitive uptake ratio 0.32-0.62).

### Where the oxidised As(III) goes (silica runs, alkalinity 1)

Root-zone (0-0.2 m) As inventories, season means, mmol m-2, control / silica:

| Goethite | As(III) | As(V) | Total As | Below 0.2 m | Net transport into column |
|---|---|---|---|---|---|
| 0.5% | 19.97 / 17.50 | 18.57 / 21.01 | 38.54 / 38.51 | 52.80 / 52.83 | 0.436 / 0.432 |
| 2% | 20.59 / 18.16 | 18.47 / 20.90 | 39.06 / 39.06 | 52.85 / 52.86 | 0.444 / 0.444 |
| 6% | 21.20 / 18.75 | 18.52 / 20.97 | 39.72 / 39.72 | 53.63 / 53.63 | 0.446 / 0.446 |

- The As(III) "lost" with silica (about 2.45 mmol m-2) reappears as As(V) in the same root zone. Total As, the deeper
  inventory and exchange across the column boundaries do not change. It is a speciation shift, not an export.
- The protective effect therefore holds only while that As(V) stays oxidised and sorbed. It can be re-reduced next
  season; the year 1 vs year 3 comparison in the sweep shows whether it persists.
- It rests on the uncalibrated As(III) oxidation constant (`as3ox`/`as3oa` 1000 /yr, O2 half-saturation 2e-5 M).

### Pre-registered correction for the 27-soil sweep (fixed 2 Nov 2026, before the sweep results)

From the silica runs only (8 soils): R_inv = 0.84 (range 0.83-0.85); log10(R_As / rule) = -0.145 (range -0.11 to -0.17).
Two corrected predictors will be scored against the basalt sweep without refitting:

- A (inventory only): R_As = 0.84 / ((1-f) + f/q)
- B (full offset): R_As = 10^-0.145 / ((1-f) + f/q) = 0.72 / ((1-f) + f/q)

Score: RMSE of log10 R_As and the share of soils whose uptake class (with T) is unchanged, against the uncorrected rule.
Basalt adds pH, Fe and redox effects the silica runs lack, so this is an out-of-sample test.

## Grain residual structure (15 treatments, `paddy/arsenic/grain_residuals.py`)

Residual = log10(observed grain ratio) - log10(competitive uptake ratio, Ki 300 uM). No global fit; descriptive only.

- Overall: mean +0.01, sd 0.17 log10 (about a factor of 1.5). No bias by growth-stage window, porewater As species
  (As(III) vs total) or grain As species (iAs vs total).
- By treatment: Si gel -0.10 (grain fell more than competition predicts; n = 4), Si-rich amendments 0.00 (n = 6), slag
  -0.05 (n = 2), diatomite +0.23 (n = 3; little Si, porewater As near detection, grain rose in two of three).
- Against porewater Si: Spearman rho -0.31 (p 0.27). The sign fits Si down-regulation of uptake at high Si, but 15
  treatments cannot resolve it.
- Grain vs straw (n = 10): grain fell less than straw in 8 of 10 (log10 grain - straw +0.13 median); slag is the
  exception (both rose together). Si cuts shoot As more than grain As, so grain is partly buffered from the uptake change
  (remobilisation or node partitioning), the place to add OsABCC1 sequestration separately.
- Across all 15, log10 grain = 0.53 x log10 predicted uptake - 0.04 (r 0.69). Grain moves about half as much as predicted
  uptake in log terms; part of that is attenuation from noisy porewater data, so treat 0.5 as a lower bound on the slope.

## 27-soil basalt sweep: first results (Stanage job 12049941, 3 Nov 2026)

Year 3, 20 t/ha basalt, all 27 soils (soils 11, 17 and 20 rerun after intermittent failures). Data:
`/mnt/project-files/arsenic-plan/sweep/as_sweep_summary_20261103.csv`; figures from `as_sweep_figures.py`.

- All 27 soils: uptake suppressed. Competitive As(III) uptake ratio (Ki 300 uM) 0.14-0.51 (Ki 600: at most 0.74); the
  simulated as3up ratio agrees (0.14-0.51). Basalt raised Si about 12.5x (T 5.8-6.0) and pH by 0.7-1.4; As(III) rose 0.95-3.4x.
- Ferrihydrite share f orders As(III) mobilisation (Spearman 0.98). The pre-registered correction A (0.84 x rule)
  predicts R_As with log RMSE 0.026 (bias +0.01) against 0.088 for the uncorrected rule; B (0.72x) over-corrects
  (RMSE 0.063). R_inv 0.80-0.88.
- Least suppression at high ferrihydrite, low goethite (3% / 0.5%: 0.38-0.51). Irrigation alkalinity has a small
  effect (alkaline water slightly less suppression, e.g. 0.33 -> 0.44 at 1.2% / 0.5%).
- To raise uptake, R_As would need to exceed T (about 5.9). Even the uncorrected rule gives at most 0.60 x T in this
  grid, so no soil here is near the boundary.
- Years 1-2 are pruned to porewater only (KEEPPH); the year 1 vs 3 comparison needs the flux-based uptake ratio.

### Year 1 vs 3 and grain iAs (3 Nov 2026)

- Year by year (`as_sweep_years.py`, from the as3up fluxes kept in pruned years): uptake ratio 0.17-0.54 in year 1,
  0.14-0.51 in year 3; it falls from year 1 to 3 in all 27 soils (median -11%). Suppression strengthens with time.
- Grain (`grain_from_uptake.py`): log10 grain = -0.037 + 0.53 log10 uptake (15 published treatments, r 0.69) gives a
  grain iAs ratio of 0.32-0.64 across the 27 soils (36-68% lower), bootstrap 5-95% range 0.21-0.76. The fitted line comes
  from Si amendments, not basalt, and the slope is attenuated by noisy porewater data; grain following uptake one-for-one
  (0.14-0.51) is the other bound.

## Next

1. Calibrate `f_site`, `Sx` and the redox constants against porewater As time series. Candidates are Seyfferth & Fendorf 2012, He 2020 and Limmer 2020.
2. Run from the settled v25 states (Stanage).
3. Calibrate `as3up` and add OsABCC1 vacuolar sequestration for grain partitioning.
4. Step 2: the Dzombak and Morel diffuse-layer model in the code, verified against PHREEQC example 8.
