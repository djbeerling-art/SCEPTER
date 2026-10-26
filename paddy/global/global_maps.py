# Maps of the global archetype results (collect.py output): per-ha season change vs no basalt in each rice cell.
# Usage: python3 paddy/global/global_maps.py <outdir with cell_ghg.csv> <out.png> [doses, default 10,20]
# Rows: dose; columns: dCH4, net realised (dCH4 + dN2O - alkalinity export - calcite), net if all released base
# cations count as CDR. Diverging scale centred on 0; red = more emissions (a penalty), blue = a net benefit.
import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, cartopy.crs as ccrs, cartopy.feature as cf

D, out = sys.argv[1], sys.argv[2]
doses = [float(x) for x in (sys.argv[3] if len(sys.argv) > 3 else '10,20').split(',')]
c = pd.read_csv(f'{D}/cell_ghg.csv')
for k in ['d_ch4', 'net', 'net_pot']: c[k] = c[k+'_tot']/c.area      # t CO2e/ha
COLS = [('d_ch4', 'ΔCH$_4$', 1.5), ('net', 'Net realised (one season)', 1.5), ('net_pot', 'Net with all released cations as CDR', 4)]
fig, axs = plt.subplots(len(doses), 3, figsize=(15, 3.3*len(doses)+0.6), dpi=200,
                        subplot_kw=dict(projection=ccrs.PlateCarree()), squeeze=False)
for i, d in enumerate(doses):
    s = c[c.dose == d]
    for j, (k, t, vm) in enumerate(COLS):
        ax = axs[i, j]; ax.set_extent([-125, 155, -35, 55], ccrs.PlateCarree())
        ax.add_feature(cf.LAND, facecolor='#f2f1ee', zorder=0); ax.add_feature(cf.COASTLINE, linewidth=0.3, edgecolor='#8a8986')
        o = s.sort_values('area')
        sc = ax.scatter(o.lon, o.lat, c=o[k], s=0.4, cmap='RdBu_r', vmin=-vm, vmax=vm, linewidths=0, transform=ccrs.PlateCarree())
        w = np.average(s[k], weights=s.area)
        ax.set_title(f'{t}, {d:.0f} t/ha  (area-weighted mean {w:+.2f})', fontsize=8.5)
        cb = fig.colorbar(sc, ax=ax, orientation='horizontal', fraction=0.05, pad=0.03, extend='both'); cb.ax.tick_params(labelsize=7)
        cb.set_label('t CO$_2$e ha$^{-1}$ season$^{-1}$ vs no basalt', fontsize=7)
fig.suptitle('Basalt on ricelands within 250 km of basaltic outcrops: one-season GHG change per rice cell (model v5, 201 archetypes)', fontsize=10)
fig.tight_layout(); fig.savefig(out)
