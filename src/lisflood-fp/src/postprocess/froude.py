#!/usr/bin/env python3
from osgeo import gdal
import math
import sys
import numpy as np

def main():
    if len(sys.argv) <= 2:
        print('Usage: froude.py <wd_in> <speed_in>')
        sys.exit(1)

    g = 9.80665
    tolh = 1e-1

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    wd_file = gdal.Open(sys.argv[1], gdal.GA_ReadOnly)
    speed_file = gdal.Open(sys.argv[2], gdal.GA_ReadOnly)
    wd_band = wd_file.GetRasterBand(1)
    speed_band = speed_file.GetRasterBand(1)

    wd_field = wd_band.ReadAsArray()
    speed_field = speed_band.ReadAsArray()

    transform = wd_file.GetGeoTransform()
    wd_nodata = wd_band.GetNoDataValue()
    speed_nodata = speed_band.GetNoDataValue()

    print("ncols", min(wd_field.shape[1], wd_field.shape[1]))
    print("nrows", min(wd_field.shape[0], wd_field.shape[0]))
    print("xllcorner", transform[0])
    print("yllcorner", transform[3]-transform[1]*wd_field.shape[0])
    print("cellsize", transform[1])
    print("NODATA_value", wd_nodata)

    for wd_row, speed_row in zip(wd_field, speed_field):
        for wd, speed in zip(wd_row, speed_row):
            if wd == wd_nodata or speed == speed_nodata:
                print(wd_nodata, end=' ')
            else:
                if wd < tolh:
                    print('0.0', end=' ')
                else:
                    print("{:.6f}".format(speed/math.sqrt(g*wd)), end=' ')
        print()

    wd_file = None
    speed_file = None

if __name__ == '__main__':
    main()
