#!/usr/bin/env python3
from osgeo import gdal
import math
import sys
import numpy as np

def main():
    if len(sys.argv) <= 2:
        print('Usage: diff.py <A_in> <B_in>', file=sys.stderr)
        sys.exit(1)

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    A_file = gdal.Open(sys.argv[1], gdal.GA_ReadOnly)
    B_file = gdal.Open(sys.argv[2], gdal.GA_ReadOnly)
    A_band = A_file.GetRasterBand(1)
    B_band = B_file.GetRasterBand(1)

    A_field = A_band.ReadAsArray()
    B_field = B_band.ReadAsArray()

    transform = A_file.GetGeoTransform()
    A_nodata = A_band.GetNoDataValue()
    B_nodata = B_band.GetNoDataValue()

    print("ncols", min(A_field.shape[1], B_field.shape[1]))
    print("nrows", min(A_field.shape[0], B_field.shape[0]))
    print("xllcorner", transform[0])
    print("yllcorner", transform[3]-transform[1]*A_field.shape[0])
    print("cellsize", transform[1])
    print("NODATA_value", A_nodata)

    for A_row, B_row in zip(A_field, B_field):
        for A, B in zip(A_row, B_row):
            if (A == A_nodata or B == B_nodata):
                print(A_nodata, end=' ')
            else:
                print("{:.20f}".format(A - B), end=' ')
        print()

    A_file = None
    B_file = None

if __name__ == '__main__':
    main()

