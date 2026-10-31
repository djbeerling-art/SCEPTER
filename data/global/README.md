# Global riceland inputs (2026-10-25)

Code: paddy/global/ on branch claude/project-thread-k5u94t (rebased on paddy/fe3-reduction).

## Files
- `cells.csv`: the run domain. 109,319 rice cells (5 arcmin), 72.4 Mha (45% of global rice area), within 250 km of
  GLiM basic volcanics. Columns: cell, lon, lat, iso3, area_irr, area_rf, ph, cec, soc, clay, temp, dist_basalt (km),
  dist_rock (km to the nearest basalt analysis), rock (data/dust_glbNN.in).
- `cells_all.csv`: all 265,876 MapSPAM rice cells, 160.7 Mha; rock = hillhouse placeholder; no basalt filter.
- `rocks.csv` (lon, lat, rock of each analysis used), `rocks_groups.csv` (median oxides and recipe of each group).
- `glim_vb.csv`: GLiM basic-volcanic grid points used for the 250 km domain.
- `arch/archetypes.csv`, `arch/cell_map.csv`: 201 archetypes (116 irrigated, 85 rainfed) from cells.csv.

## Sources and processing
- Rice area: MapSPAM 2010 v2.0 (Harvard Dataverse doi:10.7910/DVN/PRFF8V), H_TA and H_TI; rainfed = TA - TI.
  SPAM 2020 sits behind a Dataverse guestbook form, so it was not fetched.
- Soil: SoilGrids 2.0 5 km aggregates: pH(H2O), CEC (cmol/kg), SOC (wt %), clay (wt %);
  0-15 cm = (5*[0-5] + 10*[5-15])/15. Masked pixels filled from the nearest valid ring within 3 pixels.
- Temperature: WorldClim 2.1 5-arcmin tavg, mean of the warmest 4 consecutive months (no crop calendar).
- Basalt domain: GLiM (Hartmann & Moosdorf 2012), the 0.5 degree raster distributed via PANGAEA
  (dominant class per cell), class vb. Small outcrops are lost at 0.5 degree, so 72.4 Mha is a lower bound.
  The full-resolution vector GLiM would fix this if supplied.
- Rock chemistry: GEOROC (georoc.eu), subaerial whole-rock analyses with basalt rock names, downloaded in 10 degree
  tiles around the domain (258,692 rows). Kept: SiO2 45-52 and Na2O+K2O <= 5 (volatile-free TAS basalt field),
  LOI <= 3, total incl. LOI 97-101.5, not flagged altered: 13,207 analyses. 1,380 of them are the nearest analysis
  (within 250 km) for 50.1 Mha of domain rice. They were grouped by area-weighted k-means on major oxides into
  12 recipes (glb00-glb11), each = area-weighted median oxides normed by basalt_norm.py. Fe3+/Fetot from the
  analyses where most report FeO and Fe2O3, otherwise 0.15. Groups are compositional, not regional.
- Caveat: 22 Mha of domain rice has no GEOROC analysis within 250 km and takes the nearest one further away
  (area-weighted median dist_rock 182 km). EarthChem was not added.

## Area-weighted percentiles, all rice (5/25/50/75/95)
pH 5.1/5.6/6.0/6.7/7.7; CEC 11.8/16.3/18.8/22.9/31.8; SOC 0.96/1.51/2.04/3.03/6.02 %; temp 21.9/26.4/28.3/29.6/31.8 C;
clay 21.0/26.9/30.7/35.2/41.7 %.

## Rebuild
  python3 paddy/global/glim_basalt.py glim_wgs84_0point5deg.txt.asc glim_vb.csv --codes=2
  python3 paddy/global/prep_cells.py <data dir> cells_all.csv
  python3 paddy/global/rock_recipes.py rocks <georoc tile csvs> --cells=<domain cells> --n=12
  python3 paddy/global/prep_cells.py <data dir> cells.csv --glim=glim_vb.csv --rocks=rocks.csv
  python3 paddy/global/archetypes.py cells.csv arch --n=200
