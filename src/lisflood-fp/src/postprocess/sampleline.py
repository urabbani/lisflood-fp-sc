#!/usr/bin/env python3
from osgeo import gdal
import sys

def usage():
    print('Usage: sampleline.py [i|j]=n <raster_in>')
    sys.exit(1)

def main():
    if len(sys.argv) <= 2:
        usage()

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    raster = gdal.Open(sys.argv[2], gdal.GA_ReadOnly)
    band = raster.GetRasterBand(1)
    field = band.ReadAsArray()

    transform = raster.GetGeoTransform()
    blx = transform[0]
    tly = transform[3]
    dx = dy = transform[1]

    def coords(i, j):
        return blx + 0.5*dx + i*dx, tly - 0.5*dy - j*dy

    [direction, index] = sys.argv[1].split('=')
    index = int(index)

    if direction == 'i':
        for j in range(field.shape[0]):
            print(index, j, *coords(index, j), field[j, index])
    elif direction == 'j':
        for i in range(field.shape[1]):
            print(i, index, *coords(i, index), field[index, i])
    else:
        usage()

if __name__ == '__main__':
    main()
