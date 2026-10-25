# Basalt whole-rock analyses (GEOROC / EarthChem downloads) -> a small set of regional rock recipes for the global runs.
# 1. Read one or more CSVs, find latitude/longitude and oxide columns by name (GEOROC 'SIO2(WT%)', EarthChem 'SIO2'...).
# 2. Keep fresh basalts: SiO2 45-52 and Na2O+K2O <= 5 wt% (TAS basalt field, volatile-free), LOI <= --loi (default 3)
#    where reported, oxide total incl. LOI 97-101.5 wt%, valid coordinates; GEOROC whole rock and not 'altered'.
# 3. Keep samples within --dmax km (default 250) of a rice cell (--cells, a prep_cells.py output) and weight each
#    sample by the rice area for which it is the nearest analysis.
# 4. k-means (area-weighted, standardised oxides) into --n groups (default 12); each group's recipe is the
#    area-weighted median of its oxides, normed by basalt_norm.norm() into data/dust_glb<k>.in. Fe: FeO and Fe2O3
#    as reported where both are; FeOT/Fe2O3T only -> basalt_norm's Fe3+/Fetot = 0.15 convention.
# Output: <out>.csv (lon, lat, rock per sample, for prep_cells.py --rocks) and <out>_groups.csv (oxides + recipe).
# Usage: python3 rock_recipes.py <out prefix> <csv> [<csv> ...] --cells=cells_all.csv [--n=12] [--dmax=250] [--loi=3]
import sys, os, re
import numpy as np, pandas as pd
from sklearn.cluster import KMeans
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from basalt_norm import norm, write
from prep_cells import nearest, opt

OX = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'FeO', 'FeOT', 'Fe2O3T', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'LOI']
FEAT = ['SiO2', 'TiO2', 'Al2O3', 'FeOT', 'MgO', 'CaO', 'Na2O', 'K2O']

ALIAS = {**{o.lower(): o for o in OX}, 'feotot': 'FeOT', 'feototal': 'FeOT', 'fe2o3tot': 'Fe2O3T',
         'fe2o3total': 'Fe2O3T', 'latitude': 'lat', 'lat': 'lat', 'latitudemean': 'lat', 'longitude': 'lon',
         'lon': 'lon', 'long': 'lon', 'longitudemean': 'lon'}

def canon(col):
    """Column name -> oxide / lat / lon, ignoring case, units and punctuation ('SIO2(WT%)', 'Fe2O3T', 'LATITUDE MEAN')."""
    return ALIAS.get(re.sub(r'\(.*?\)|wt%|[^a-z0-9]', '', col.lower()))

def read(path):
    """One download -> lat, lon and oxide columns. GEOROC: lat/lon = mean of MIN/MAX, whole-rock analyses only,
    samples flagged altered dropped."""
    d = pd.read_csv(path, low_memory=False, encoding='latin-1', index_col=False, on_bad_lines='skip')   # GEOROC appends references
    up = {c.upper(): c for c in d.columns}
    if 'TYPE OF MATERIAL' in up: d = d[d[up['TYPE OF MATERIAL']].astype(str).str.upper().str.contains('WHOLE ROCK')]
    if 'ALTERATION' in up: d = d[~d[up['ALTERATION']].astype(str).str.upper().str.contains('ALTERED')]
    if 'LATITUDE (MIN.)' in up:
        for k in ('LATITUDE', 'LONGITUDE'):
            mm = d[[up[f'{k} (MIN.)'], up.get(f'{k} (MAX.)', up[f'{k} (MIN.)'])]].apply(pd.to_numeric, errors='coerce')
            d = d.assign(**{k.lower(): mm.mean(axis=1)})
    ren = {}
    for c in d.columns:
        k = canon(c)
        if k and k not in ren.values(): ren[c] = k
    d = d[list(ren)].rename(columns=ren).apply(pd.to_numeric, errors='coerce')
    for o in OX + ['lat', 'lon']:
        if o not in d: d[o] = np.nan
    return d

def clean(d, loi):
    d = d.dropna(subset=['lat', 'lon', 'SiO2', 'MgO', 'CaO', 'Al2O3'])
    d = d[d.lat.between(-90, 90) & d.lon.between(-180, 180)].copy()
    feot = d.FeOT.fillna(d.Fe2O3T*0.8998)
    both = d.FeO.notna() & d.Fe2O3.notna()
    feot[both] = d.FeO[both] + 0.8998*d.Fe2O3[both]
    feot = feot.fillna(d.FeO)
    d['FeOT'] = feot
    d.loc[~both, ['FeO', 'Fe2O3']] = np.nan
    d = d.dropna(subset=['FeOT'])
    ox = ['SiO2', 'TiO2', 'Al2O3', 'FeOT', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5']
    tot = d[ox].fillna(0).sum(axis=1)
    d = d[(tot + d.LOI.clip(lower=0).fillna(0)).between(97, 101.5) & (d.LOI.isna() | (d.LOI <= loi))]
    vf = d[ox].fillna(0).div(d[ox].fillna(0).sum(axis=1), axis=0)*100       # volatile-free for TAS
    keep = vf.SiO2.between(45, 52) & (vf.Na2O + vf.K2O <= 5)
    return d[keep].drop_duplicates(subset=['lat', 'lon'] + ox)

def wmedian(x, w):
    m = ~np.isnan(x)
    if not m.any(): return np.nan
    o = np.argsort(x[m]); cw = np.cumsum(w[m][o])
    return x[m][o][np.searchsorted(cw, cw[-1]/2)]

if __name__ == '__main__':
    out = sys.argv[1]
    files = [a for a in sys.argv[2:] if not a.startswith('--')]
    n, dmax, loi = int(opt('n', 12)), float(opt('dmax', 250)), float(opt('loi', 3))
    s = clean(pd.concat([read(f) for f in files], ignore_index=True), loi)
    print(f'{len(s)} fresh basalt analyses with coordinates')
    c = pd.read_csv(opt('cells', 'cells_all.csv'))
    dist, i = nearest(c.lon.values, c.lat.values, s.lon.values, s.lat.values)
    area = np.bincount(i[dist <= dmax], weights=(c.area_irr+c.area_rf).values[dist <= dmax], minlength=len(s))
    s = s[area > 0].assign(area=area[area > 0]).reset_index(drop=True)
    print(f'{len(s)} analyses are the nearest within {dmax:.0f} km to {s.area.sum()/1e6:.1f} Mha of rice')
    X = s[FEAT].fillna(s[FEAT].median()); X = ((X - X.mean())/X.std().replace(0, 1)).values
    k = min(n, len(s))
    s['grp'] = KMeans(k, n_init=10, random_state=0).fit(X, sample_weight=s.area.values).labels_
    rows = []
    for g, m in s.groupby('grp'):
        w = m.area.values
        ox = {o: wmedian(m[o].values, w) for o in ['SiO2', 'TiO2', 'Al2O3', 'FeOT', 'MnO', 'MgO', 'CaO', 'Na2O',
                                                   'K2O', 'P2O5', 'FeO', 'Fe2O3']}
        name = f'glb{g:02d}'
        oxn = {o: v for o, v in ox.items() if not np.isnan(v) and o not in ('FeO', 'Fe2O3')}
        if (m.FeO.notna() & m.Fe2O3.notna()).mean() > 0.5:      # Fe split from the analyses that report it
            r3 = wmedian((0.8998*m.Fe2O3/(m.FeO + 0.8998*m.Fe2O3)).values, w)
            oxn['Fe2O3'] = r3*ox['FeOT']/0.8998; oxn['FeO'] = (1-r3)*ox['FeOT']; oxn.pop('FeOT')
        r = norm(oxn)
        write(name, r, f'area-weighted median of {len(m)} basalt analyses (paddy/global/rock_recipes.py)')
        rows.append(dict(rock=name, nsample=len(m), area=w.sum(), lon=np.average(m.lon, weights=w),
                         lat=np.average(m.lat, weights=w), **{o: round(v, 2) for o, v in ox.items()},
                         recipe=' '.join(f'{k}:{v:.3f}' for k, v in sorted(r.items(), key=lambda kv: -kv[1]))))
    s['rock'] = [f'glb{g:02d}' for g in s.grp]
    s[['lon', 'lat', 'rock']].to_csv(f'{out}.csv', index=False)
    pd.DataFrame(rows).to_csv(f'{out}_groups.csv', index=False)
    print(pd.DataFrame(rows)[['rock', 'nsample', 'area', 'lat', 'lon', 'SiO2', 'MgO', 'FeOT', 'recipe']].to_string())
