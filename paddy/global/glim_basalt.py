# GLiM (Hartmann & Moosdorf 2012, doi:10.1594/PANGAEA.788537) -> points on a 5-arcmin grid covering the chosen
# lithology classes, for prep_cells.py --glim (distance from each rice cell to the nearest outcrop).
#   vector input (LiMW_GIS 2015.gdb or any OGR file with the first-level class in field 'xx'): polygons rasterised
#   with all_touched, so small outcrops and dykes are kept; raster input (GLiM ASCII/GeoTIFF of class codes):
#   pass the codes with --codes, since their numbering depends on the release.
# Default class: vb (basic volcanics). --classes=vb,pb also takes basic plutonics (gabbro).
# Usage: python3 glim_basalt.py <glim file> <out csv> [--classes=vb] [--codes=..] [--res=0.0833333]
import sys
import numpy as np, pandas as pd

def opt(name, default):
    return next((a.split('=', 1)[1] for a in sys.argv if a.startswith(f'--{name}=')), default)

def grid_points(mask, transform):
    i, j = np.nonzero(mask)
    lon, lat = transform * (j + 0.5, i + 0.5)
    return pd.DataFrame({'lon': np.round(lon, 4), 'lat': np.round(lat, 4)})

def from_vector(path, classes, res):
    import geopandas as gpd, rasterio.features
    from rasterio.transform import from_origin
    g = gpd.read_file(path, columns=['xx'])
    g = g[g.xx.str.strip().str.lower().isin(classes)].to_crs('EPSG:4326')
    print(f'{len(g)} polygons of class {",".join(classes)}')
    tr = from_origin(-180, 90, res, res)
    shape = (int(round(180/res)), int(round(360/res)))
    m = rasterio.features.rasterize(((geom, 1) for geom in g.geometry if geom is not None), out_shape=shape,
                                    transform=tr, all_touched=True, dtype='uint8')
    return grid_points(m > 0, tr)

def from_raster(path, codes):
    import rasterio
    with rasterio.open(path) as r:
        a = r.read(1)
        if r.crs is not None and r.crs.to_epsg() not in (None, 4326):
            sys.exit('raster GLiM must be in geographic lon/lat')
        return grid_points(np.isin(a, codes), r.transform)

if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    codes = opt('codes', None)
    if codes:
        p = from_raster(src, [int(x) for x in codes.split(',')])
    else:
        p = from_vector(src, opt('classes', 'vb').lower().split(','), float(opt('res', 1/12)))
    p.to_csv(out, index=False)
    print(f'wrote {out}: {len(p)} outcrop grid points')
