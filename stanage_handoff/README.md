# Handoff: global riceland SCEPTER runs on Stanage (31 Oct 2026)

Goal (David, 31 Oct 2026): run the global riceland enhanced-weathering simulations on the University of Sheffield
Stanage HPC, using the paddy SCEPTER model as frozen for the GMD paper (v22). These runs are for the follow-on
global paper (Nature Geoscience target), not the GMD paper, which has no global runs.

## 1. Access (working as of 31 Oct 07:03 UTC)
- A Remote Control session on David's university PC (UOSD4A2CD73110B, Windows 11, folder
  `C:\Users\bo1djb\Documents\Claude Oct 2026`) can run `ssh bo1djb@stanage.shef.ac.uk` with no password or MFA,
  using the key `C:\Users\bo1djb\.ssh\id_ed25519` (no passphrase). It needs David's FortiClient VPN to be connected.
- Login node: logina.stanage.hpc.shef.ac.uk. Slurm account `default`; partition `sheffield` (4-day limit) for
  production, `interactive` (8 h) for tests. Storage: 50 GB backed-up home; scratch `/mnt/parscratch/users/bo1djb`
  (linked as `~/fastdata`, not backed up) for run output. Modules include GCC, Python 3.11-3.13, netCDF, GDAL;
  no conda (OpenBLAS module not yet checked). Survey done 31 Oct from the PC session; a parallel starter kit
  (HANDOFF.md, PowerShell copy/test scripts, generic Slurm driver) is on the PC in
  `C:\Users\bo1djb\Documents\Claude Oct 2026\riceland-stanage\`. Its open questions (code, rice map, inputs,
  scenarios) are answered by this package; reuse its PowerShell transfer scripts.
- The cloud container cannot reach Stanage directly. Code moves via GitHub (`git clone` on Stanage); results move back
  via the PC session (scp to the PC) or by pushing small summary files to GitHub.

## 2. Model version
**DECIDED (David, 31 Oct 07:24 UTC): v23 replaces v22.** Executable `scepter_paddy_v23` from commit 267a856
(v22 + carbon-closure fix, 9 lines in scepter.f90); options and `krxn_calibrated_v22.in` unchanged; `v22_options.env`
now sets `EXE=scepter_paddy_v23`; `slurm/build.sh` builds it. Recalibrate only if the GMD thread reports a material shift.
The v22 notes below are kept for provenance.

- Repository: github.com/djbeerling-art/SCEPTER, branch `claude/project-thread-lst1mi` (the global-run thread; = `claude/project-thread-pip6bj-pan` + archetypes_v22, pilot list, collect_chains.py; this package is in `stanage_handoff/`).
- Executable: `scepter_paddy_v22` (source = commit 162913a, unchanged up to 6e71f12), built from `makefile_paddy`:
  `make --file=makefile_paddy PROGRAM=scepter_paddy_v22 OBJS=scepter_paddy_v22.o` (gfortran + OpenBLAS; see
  `slurm/build.sh`). Python: numpy, scipy, pandas, scikit-learn, SALib.
- Calibration constants: `paddy/krxn_calibrated_v22.in` (= the cloud `k_cal.in`; g2ch4 vmax 1860).
- Exact options: `v22_options.env` (SOIL for spin-ups, PROD for season/fallow chains, P21 = 0.5 m column with
  plough pan, TILL=0.2). Si feedback K_Si 2 mM, w 0.2.
- PENDING: a carbon-closure fix (v23) is being tested in the GMD thread. Upstream SCEPTER OM decay ("g2") sends 4% of
  its C to oxalate/acetate, which the paddy runs do not track, so 0.35% of gross C was lost. v23 returns that C to
  CO2. Do not start production global runs until David decides whether v23 replaces v22; the setup (build,
  test archetypes, Slurm plumbing) can go ahead on v22 now, and switching is one variable (`EXE`).

## 3. Inputs (in /mnt/project-files/riceland-scepter/global/, see its README.md)
- `cells.csv`: 109,319 rice cells (5 arcmin), 72.4 Mha within 250 km of GLiM basic volcanics (lower bound; 0.5 deg
  raster). Columns include area_irr, area_rf, ph, cec, soc, clay, temp, rock (glb00-glb11 recipes, data/dust_glbNN.in).
- `cells_all.csv`: all 265,876 MapSPAM 2010 cells, 160.7 Mha (rock = Hillhouse placeholder).
- `arch/archetypes.csv` + `cell_map.csv`: 201 archetypes (116 irrigated, 85 rainfed) from `paddy/global/archetypes.py`.
- Rock recipes `data/dust_glb00.in` ... `glb11.in` are in the repository.

### Must fix before production: the soil-conditioning rule in archetypes.py is out of date
`archetypes.py` maps target soil pH to irrigation alkalinity (`--irralk`) and soil calcite (`--soilcc`) with a lookup
built on the v5 15-cell column (Oct 25). The v22 runs (GMD archetypes, runR.sh, 29 Oct) used a different site rule:
- dry-soil pH >= 6.8: `--soilcc=0.005 --irralk=10`;
- otherwise `--irralk = interp(pH; [5.0, 6.0, 6.2, 6.5] -> [0, 1, 4, 10])`.
DONE (31 Oct): `stanage_handoff/arch/archetypes_v22.csv` = the same 201 archetypes (ids and membership unchanged,
`cell_map.csv` copied alongside) with args recomputed by `python3 paddy/global/archetypes.py <old csv> --reargs <new csv>`
(rule now in `ph_to_soil`; 45 archetypes, 21% of area, are calcareous). It reproduces the GMD args to rounding
(irr061 pH 6.045 gives irralk 1.67 vs 1.6 used in the GMD run, rf189 6.8 exactly). Still to check: settled porewater
pH on a few archetypes (pilot). Rainfed archetypes get `--awd` (the 14-day pre-harvest drain applies to both regimes).

## 4. Run protocol (identical to the GMD v22 archetype runs)
Per archetype (`run_arch_v22.sh`):
1. 3-yr drained spin-up (`run_flooded_column.py ... --fallow --basalt=0 $SOIL $P21 <args>`).
2. 7-yr settle without rock (`paddy/multiyear/run_chain.sh ... 7 0`), keeping only the year-7 fallow as the restart.
3. Per dose (0, 10, 20 t/ha/yr), copy the settled state and run an N-year chain (113-d season with basalt at its
   start + 252-d drained fallow). Treatment minus the dose-0 chain (year-matched) gives the effect.
4. Budgets: `paddy/ghg_budget.py` per year (season + fallow), with PDA (= ALK - 2Fe(II) - 2NH4+) export and calcite
   capped at cations not exported; GWP100 (CH4 27, N2O 273). For the GMD archetypes this was `fig8_arch_v22.py` /
   `figcdr22`; `paddy/global/collect.py` is the old one-season version; for chains use `paddy/global/collect_chains.py`
   (same budget as figcdr22, PDA by default; writes arch_years/arch_cum/cell_cum/nation_cum CSVs, small enough to push).

Cost (cloud container, 1 core): spin-up ~3 min; settle ~35 min; a 10-yr chain 1.5 h (control) to ~3.5 h (20 t/ha;
basalt years slow to ~25 min after year 5). About 8-9 core-hours per archetype for 0/10/20 t over 10 yr, so 201
archetypes is ~1,800 core-hours: one array job (`slurm/arch_array.sbatch`, one task per archetype, 24 h wall).
A one-season-per-dose run (no multi-year) is ~1 core-hour per archetype.
Disk: ~30 MB per pruned year; ~1 GB per archetype for 3 doses x 10 yr, ~200 GB total. Keep `prof/` only for the
restart states (run_chain.sh already prunes).

## 5. Decisions for David before production
1. v22 or v23 (carbon-closure fix; see section 2).
2. Length: 10-yr chains (GMD standard) or one season per dose. David earlier preferred the multi-year route
   ("global 10-yr ~60 h on 4 cores" before Stanage; trivial on Stanage).
3. Doses: 0/10/20 t/ha/yr (GMD), optionally 30.
4. Si feedback: w 0.2 central; w 0 (feedback off) as the structural sensitivity, since the Morris and Fig 6b results
   show w sets the sign and size of the CH4 response. Recommend a w 0 global set alongside.
5. Domain: cells.csv (72.4 Mha near basalt) for production; cells_all.csv only for a context map.
6. Known input caveats: soil reducible Fe(III) is a clay proxy (0.012*clay/30, clipped 0.004-0.03), and soil
   ferrihydrite is the top control on absolute net GHG in the v22 Morris screen, so run the clip bounds as a
   sensitivity; one rice season per year (multiple cropping not counted); GLiM raster domain is a lower bound;
   22 Mha uses a rock analysis > 250 km away.

## 6. Suggested first steps in the new thread
1. On Stanage: clone the repo, run `slurm/build.sh` (confirm module names), and run one archetype for 1 year at
   dose 0 and 20 to check the build reproduces the cloud numbers (the GMD irr061 archetype: `rf/irr` args in
   `archR.txt` below).
2. (done) archetypes regenerated with the v22 soil rule; collect_chains.py written (untested until pilot output exists).
3. Submit a 10-archetype pilot array (`sbatch --array=1-10 stanage_handoff/slurm/arch_array.sbatch stanage_handoff/arch/pilot10.csv`:
   the six GMD archetypes plus the four largest others, 7.1 Mha), check timing and disk, then the full array once David picks v22/v23.

GMD archetype args used for the v22 runs (for the step-1 check; irr061 year-1 control and 20 t results are in the
GMD thread):
```
irr061|--rock=glb04 --irralk=1.6 --cec=16.9 --temp=28.7 --om=506 --fe3=0.0112
rf189|--awd --rock=glb09 --irralk=6.8 --cec=17.5 --temp=29.8 --om=533 --fe3=0.0125
```

## Files here
- `README.md` (this file), `v22_options.env`, `krxn_calibrated_v22.in`
- `run_arch_v22.sh`: one archetype, full protocol, resumable via its status file
- `slurm/build.sh`, `slurm/arch_array.sbatch`: build and array-job templates (module names unverified)
