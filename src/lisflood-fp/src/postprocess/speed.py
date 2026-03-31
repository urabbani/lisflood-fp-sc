#!/usr/bin/env python3
from osgeo import gdal
import math
import sys
import numpy as np

def main():
    if len(sys.argv) <= 3:
        print('Usage: speed.py <Vx_in> <Vy_in> <speed_out>')
        sys.exit(1)

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    Vx_file = gdal.Open(sys.argv[1], gdal.GA_ReadOnly)
    Vy_file = gdal.Open(sys.argv[2], gdal.GA_ReadOnly)
    Vx_band = Vx_file.GetRasterBand(1)
    Vy_band = Vy_file.GetRasterBand(1)

    Vx_field = Vx_band.ReadAsArray()
    Vy_field = Vy_band.ReadAsArray()

    transform = Vx_file.GetGeoTransform()
    nodata = Vx_band.GetNoDataValue()

    with open(sys.argv[3], 'w') as out:
        print("ncols", min(Vx_field.shape[1], Vy_field.shape[1]), file=out)
        print("nrows", min(Vx_field.shape[0], Vy_field.shape[0]), file=out)
        print("xllcorner", transform[0], file=out)
        print("yllcorner", transform[3]-transform[1]*Vx_field.shape[0], file=out)
        print("cellsize", transform[1], file=out)
        print("NODATA_value", nodata, file=out)

        for Vx_row, Vy_row in zip(Vx_field, Vy_field):
            for Vx, Vy in zip(Vx_row, Vy_row):
                if (Vx == nodata):
                    print(nodata, end=' ', file=out)
                else:
                    print("{:.20f}".format(math.sqrt(Vx**2 + Vy**2)), end=' ', file=out)
            print(file=out)

    Vx_file = None
    Vy_file = None

if __name__ == '__main__':
    main()
