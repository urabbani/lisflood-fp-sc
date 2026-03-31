#!/usr/bin/env python3
from osgeo import gdal
import sys

def main():
    if len(sys.argv) <= 2:
        print('Usage: mask.py <rasterfile> <demfile>', file=sys.stderr)
        sys.exit(1)

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')

    raster = gdal.Open(sys.argv[1])
    transform = raster.GetGeoTransform()
    raster_band = raster.GetRasterBand(1)
    raster_nodata = raster_band.GetNoDataValue()
    raster_field = raster_band.ReadAsArray()

    dem = gdal.Open(sys.argv[2])
    dem_band = dem.GetRasterBand(1)
    dem_nodata = dem_band.GetNoDataValue()
    dem_field = dem_band.ReadAsArray()

    print("ncols", raster_field.shape[1])
    print("nrows", raster_field.shape[0])
    print("xllcorner", transform[0])
    print("yllcorner", transform[3]-transform[1]*raster_field.shape[0])
    print("cellsize", transform[1])
    print("NODATA_value", raster_nodata)

    for raster_row, dem_row in zip(raster_field, dem_field):
        for raster_val, dem_val in zip(raster_row, dem_row):
            if (dem_val == dem_nodata):
                print(raster_nodata, end=' ')
            else:
                print("{:.10f}".format(raster_val), end=' ')
        print()

if __name__ == '__main__':
    main()
