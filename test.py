import pluto
import numpy as np
import pyproj
from scipy import interpolate
import matplotlib.pyplot as plt
import wrfutils as wu
from scipy.interpolate import griddata


mypluto = pluto.Pluto()

data = mypluto.computeUCPs()

keys = data[0].keys()

source_proj = pyproj.Proj(init='EPSG:2263')
target_proj = pyproj.Proj(init='EPSG:4326')
conv = .3048  # foot to meter conversion factor

# Use pyproj to transform from source projection (I mean, WTF is
# New York/Long Island State Plane projection anyways? Ugh).
coords = [pyproj.transform(source_proj,
                           target_proj,
                           conv*float(x),
                           conv*float(y))
          for x, y in keys]
lon, lat = zip(*coords)

# To generate the regular grid, we use the min and max values from
# lon/lat
delta = 0.001
longrid = np.arange(min(lon), max(lon) + delta, delta)
latgrid = np.arange(min(lat), max(lat) + delta, delta)

newx, newy = np.meshgrid(longrid, latgrid)

# newfrac = griddata(mypoints, bfrac, (newx, newy), method='nearest')
# # newfrac[newfrac < 1] = np.nan
#
# fig, ax = plt.subplots()
# mp = wu.forecastMap(ax)
#
# im = mp.contourf(newx, newy, newfrac, np.arange(0, 1.01, .1),
#                  latlon=True,
#                  cmap='Spectral_r',
#                  extend='both')
# mp.readshapefile('nycd_geocoords', 'nycd')
# mp.colorbar(im)
# fig.savefig('interp-test-bfrac.png', dpi=200, bbox_inches='tight')
