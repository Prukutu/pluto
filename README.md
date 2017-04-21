# uwrf-pluto
Tools to extract urban canopy parameters from NYC's PLUTO data and interpolate to regular grids.

## Installing and using
### External dependencies

1. NumPy - Pluto methods use numpy arrays internally. Some methods also return numpy arrays.
2. SciPy - Pluto uses the interpolate module from SciPy to grid from scattered points to a regular grid
3. PyProj - Needed to convert from the native Pluto projection to Lat/Lon coordinates
