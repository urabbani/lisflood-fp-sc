#!/usr/bin/env python3
from osgeo import gdal
import math
import sys
import numpy as np

def main():
    if len(sys.argv) <= 1:
        print('Usage: downsample.py <in>', file=sys.stderr)
        sys.exit(1)

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    raster = gdal.Open(sys.argv[1], gdal.GA_ReadOnly)
    band = raster.GetRasterBand(1)
    field = band.ReadAsArray()

    transform = raster.GetGeoTransform()
    nodata = band.GetNoDataValue()
    
    raster = None

    xsz = field.shape[1] // 2
    ysz = field.shape[0] // 2

    print("ncols", xsz)
    print("nrows", ysz)
    print("xllcorner", transform[0])
    print("yllcorner", transform[3]-transform[1]*ysz*2)
    print("cellsize", transform[1]*2)
    print("NODATA_value", nodata)

    for j in range(0, ysz):
        for i in range(0, xsz):
            a = field[j*2, i*2]
            b = field[j*2, i*2+1]
            c = field[j*2+1, i*2]
            d = field[j*2+1, i*2+1]

            data = [e for e in [a, b, c, d] if abs(e-nodata) >= 1e-12]
            if len(data) > 0:
                print(np.mean(data), end=' ')
            else:
                print(nodata, end=' ')
        print()

if __name__ == '__main__':
    main()
