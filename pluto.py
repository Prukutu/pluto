import numpy as np
from os import listdir
import pyproj
import csv
from scipy.interpolate import griddata

# A module to deal with PLUTO data. It should cover projecting to geographical
# coordinates, calculating the variables needed for WRF, and interpolating into
# a regular grid as class methods.
# TODO: Incorporate parks and other areas without buildings.

class Pluto:

    def __init__(self, workdir='./'):
        """ Initialize the Pluto class with working directory. The working
            directory should contain all the raw PLUTO files. Sadly, the don't
            come as a single file for NYC. One file per borough.
        """
        self.workdir = workdir

        self.boros = ('SI', 'MN', 'BK', 'BX', 'QN')
        self.filenames = {boro: '.'.join([boro, 'csv']) for boro in self.boros}

#        We should check that the files are in workdir.
        for fname in self.filenames.values():
            assertmsg = 'Pluto borough file not in workdir!'
            assert fname in listdir(self.workdir), assertmsg

    def getRawData(self):
        """ Get the coordinates of each lot in the borough files.
            Returns 2 1-D arrays containing the lat/lon of the lots.
        """

        # Utility function to grab fields from the boro files.
        def getFields(borofile,
                      xi=-11,
                      yi=-10,
                      numflrsindx=44,
                      grgareaindx=38,
                      bldgareaindx=33,
                      lotareaindx=32,
                      bldfrontindx=49,
                      luindx=28):
            """ Get raw data from the PLUTO csv files. For now, I'm only
                grabbing the data I'm interested in for computing UCPs
            """
            with open(borofile, 'r') as f:

                # Read file lines, exclude the header.
                reader = csv.reader(f)
                splitlines = [row for row in reader][1:]

                # Return list where each item is a tuple of (Xcoord, Ycoord)
                numfloors = {(l[xi], l[yi]): l[numflrsindx]
                             for l in splitlines}
                garageArea = {(l[xi], l[yi]): l[grgareaindx]
                              for l in splitlines}
                bldgArea = {(l[xi], l[yi]): l[bldgareaindx]
                            for l in splitlines}
                lotArea = {(l[xi], l[yi]): l[lotareaindx]
                           for l in splitlines}
                bldgFront = {(l[xi], l[yi]): l[bldfrontindx]
                             for l in splitlines}

            return (numfloors, garageArea, bldgArea, lotArea, bldgFront)

        # Get the data for each borough. The reason it's done this way is
        # because PLUTO comes in separate files per borough.
        data = {boro: getFields(self.filenames[boro])
                for boro in self.boros}

        # Separate into dictionaries and flatten to remove borough
        # dependency of the data. We end up with large dicts for the entire
        # city, as opposed to a dictionary of dictionaries.
        numfloors = {key: data[bigkey][0][key] for bigkey in data.keys()
                     for key in data[bigkey][0].keys()}
        garageArea = {key: data[bigkey][1][key] for bigkey in data.keys()
                      for key in data[bigkey][1].keys()}
        bldgArea = {key: data[bigkey][2][key] for bigkey in data.keys()
                    for key in data[bigkey][2].keys()}
        lotArea = {key: data[bigkey][3][key] for bigkey in data.keys()
                   for key in data[bigkey][3].keys()}
        bldgFront = {key: data[bigkey][4][key] for bigkey in data.keys()
                     for key in data[bigkey][4].keys()}

        # delete key/value pair with empty xcoord/ycoord field.
        for urbanparam in (numfloors,
                           garageArea,
                           bldgArea,
                           lotArea,
                           bldgFront):
            del(urbanparam[('', '')])

        return (numfloors, garageArea, bldgArea, lotArea, bldgFront)
        # return data

    def computeUCPs(self):
        """ Compute the UCP's needed to run BEPBEM. We use Gutierrez's
            methodology. We will compute the following UCP's:

            urban area fraction
            building height (and height histograms?)
            building surface to height ratio
        """

        # Define functions to extract the UCPs from the PLUTO raw data.

        def getBuildAreaFraction(lotArea, buildArea, numFloors, garageArea):

            """ According to Burian et al. (2008) we define the plan area
                fraction as:

                    planAreaFrac = planArea/lotArea

                where A_p is the total plan area of buildings. We compute A_p
                from PLUTO fields as:

                    A_p = (buildArea - garageArea)/numFloors
            """

            try:

                planArea = (float(buildArea)
                            - float(garageArea))/float(numFloors)

                if float(lotArea) >= planArea:
                    planAreaFrac = planArea/float(lotArea)
                else:
                    planAreaFrac = 1.0

            except ValueError:
                # print 'At least one pluto value unavailable'
                planAreaFrac = 0
            except ZeroDivisionError:
                # print 'No floor data available!'
                planAreaFrac = 0

            return planAreaFrac

        def getBuildHeight(numfloors, floorheight=5):
            """ Compute building height from the number of floors in a
                building. Floor height can be changed (default 5 m)
            """

            try:
                buildheight = float(numfloors)*floorheight
            except ValueError:
                # Just some error catching in case of missing values, usually
                # as empty strings.
                # print 'No numfloors available'
                buildheight = 0

            return buildheight

        def getBuildSurfPlanRatio(buildhgt, buildfront, buildfrac):
            """ Compute the building surface area to plan area ratio:
                    surfplanratio = (2*buildheight/buildfront + 1)*buildfrac

                INPUT:
                buildhgt: building height in meters
                buildfront: building frontage in meters
                buildfrac: building plan area fraction
            """

            try:
                buildwidth = float(buildfront)*.3048
                surfplanratio = (2*buildhgt/float(buildwidth) + 1)*buildfrac

            except ValueError:
                surfplanratio = 0
            except ZeroDivisionError:
                surfplanratio = 0

            return surfplanratio

        data = self.getRawData()
        mykeys = data[0].keys()

        # Actually compute the UCPs.
        buildAreaFrac = {key: getBuildAreaFraction(data[3][key],
                                                   data[2][key],
                                                   data[0][key],
                                                   data[1][key])
                         for key in mykeys}
        buildheight = {key: getBuildHeight(data[0][key]) for key in mykeys}

        buildsurfplanratio = {key: getBuildSurfPlanRatio(buildheight[key],
                                                         data[4][key],
                                                         buildAreaFrac[key])
                              for key in mykeys}

        return (buildAreaFrac, buildheight, buildsurfplanratio)

    def interpolateUCP(self, newgrid=None, delta=0.001, method='nearest'):

        """ This method interpolates the scattered building data from PLUTO
            into a regular grid at grid spacing coorddelta (milli-degree
            default). We also convert from the PLUTO projection (EPSG: 2263, m)
            to a geographical coordinate system (EPSG: 4326, lat/lon).

            INPUT:
            newgrid: Tuple containing n x m shape arrays for lon and lat. If
             None, computes the target grid using the mix/ma coordinates of the
            source PLUTO data.
        """

        # Get the UCPs as dictionaries with coordinates as the key
        bfrac, bhgt, bsurf = self.computeUCPs()

        # All keys are present
        keys = bfrac.keys()  # list of (Xcoord, Ycoord) tuples

        # Define the source and target projections using pyproj
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
        if newgrid is None:
            longrid = np.arange(min(lon), max(lon) + delta, delta)
            latgrid = np.arange(min(lat), max(lat) + delta, delta)

            newx, newy = np.meshgrid(longrid, latgrid)
        else:
            newx, newy = newgrid
        sourcecoords = np.array((lon, lat)).transpose()
        # Now generate the gridded UCPs using griddata
        ucp = [griddata(sourcecoords, val.values(), (newx, newy), method=method)
               for val in (bfrac, bhgt, bsurf)]

        return newx, newy, ucp
