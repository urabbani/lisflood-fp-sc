#!/usr/bin/env python3
from osgeo import gdal
import math
import sys

def main():
    if len(sys.argv) <= 1:
        print('Usage: stats.py <rasterfile>', file=sys.stderr)
        sys.exit(1)

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    raster = gdal.Open(sys.argv[1])
    band = raster.GetRasterBand(1)
    nodata = band.GetNoDataValue()

    allow_approximation = 0
    band.ComputeStatistics(allow_approximation)
    minimum = band.GetMinimum()
    maximum = band.GetMaximum()

    mae = 0.0
    rmse = 0.0

    data_cells = 0
    nodata_cells = 0
    for row in band.ReadAsArray():
        for val in row:
            if (val == nodata):
                nodata_cells += 1
            else:
                data_cells += 1
                mae += abs(val)
                rmse += val*val

    mae /= data_cells
    rmse = math.sqrt(rmse/data_cells)

    print('# min max l_infty mae rmse total_cells nodata_cells')
    print("{:.20f} {:.20f} {:.20f}".format(minimum, maximum, max(abs(minimum), abs(maximum))), end=' ')
    print("{:.20f} {:.20f}".format(mae, rmse), end=' ')
    print(raster.RasterXSize*raster.RasterYSize, nodata_cells)

if __name__ == '__main__':
    main()
