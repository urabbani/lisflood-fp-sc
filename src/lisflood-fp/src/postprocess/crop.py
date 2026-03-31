#!/usr/bin/env python3
from osgeo import gdal
import math
import sys
import numpy as np

def main():
    if len(sys.argv) <= 5:
        print('Usage: crop.py <in> <crop_west> <crop_east> <crop_north> <crop_south>')
        sys.exit(1)

    crop_w = int(sys.argv[2])
    crop_e = int(sys.argv[3])
    crop_n = int(sys.argv[4])
    crop_s = int(sys.argv[5])

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    Vx_file = gdal.Open(sys.argv[1], gdal.GA_ReadOnly)
    Vx_band = Vx_file.GetRasterBand(1)

    Vx_field = Vx_band.ReadAsArray()

    transform = Vx_file.GetGeoTransform()
    delta = transform[1]
    nodata = Vx_band.GetNoDataValue()

    print("ncols", Vx_field.shape[1] - crop_w - crop_e)
    print("nrows", Vx_field.shape[0] - crop_n - crop_s)
    print("xllcorner", transform[0] + crop_w*delta)
    print("yllcorner", transform[3]-delta*Vx_field.shape[0] + crop_s*delta)
    print("cellsize", delta)
    print("NODATA_value", nodata)

    for j, row in enumerate(Vx_field):
        if j < crop_n or j >= Vx_field.shape[0] - crop_s:
            continue
        for i, val in enumerate(row):
            if i >= crop_w and i < Vx_field.shape[1] - crop_e:
                if (val == nodata):
                    print(nodata, end=' ')
                else:
                    print("{:.20f}".format(val), end=' ')
        print()

    Vx_file = None

if __name__ == '__main__':
    main()
