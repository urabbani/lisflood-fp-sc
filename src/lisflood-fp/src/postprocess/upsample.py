#!/usr/bin/env python3
from osgeo import gdal
import math
import sys
import numpy as np

def main():
    if len(sys.argv) <= 2:
        print('Usage: upsample.py magnitude <in>', file=sys.stderr)
        sys.exit(1)

    upsample = int(sys.argv[1])

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    raster = gdal.Open(sys.argv[2], gdal.GA_ReadOnly)
    band = raster.GetRasterBand(1)
    field = band.ReadAsArray()

    transform = raster.GetGeoTransform()
    nodata = band.GetNoDataValue()
    
    raster = None

    xsz = field.shape[1] * upsample
    ysz = field.shape[0] * upsample

    print("ncols", xsz)
    print("nrows", ysz)
    print("xllcorner", transform[0])
    print("yllcorner", transform[3]-transform[1]*field.shape[0])
    print("cellsize", transform[1]/float(upsample))
    print("NODATA_value", nodata)

    for j in range(ysz):
        for i in range(xsz):
            print(field[j//upsample][i//upsample], end=' ')
        print()

if __name__ == '__main__':
    main()
