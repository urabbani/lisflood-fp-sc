#!/usr/bin/env python3
from osgeo import gdal
import sys

def main():
    if len(sys.argv) <= 3:
        print('Usage: point.py <file> <x> <y>')
        sys.exit(1)

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    f = gdal.Open(sys.argv[1], gdal.GA_ReadOnly)
    band = f.GetRasterBand(1)
    field = band.ReadAsArray()

    x = float(sys.argv[2])
    y = float(sys.argv[3])

    transform = f.GetGeoTransform()
    ysz = field.shape[0]
    xsz = field.shape[1]
    dx = dy = transform[1]

    blx = transform[0]
    tlx = blx + dx*xsz
    tly = transform[3]
    bly = tly - dy*ysz

    if x < blx or x >= tlx or y <= bly or y > tly:
        print(f"({x}, {y}) is outside the domain ({blx}, {bly}) -- ({tlx}, {tly})")
        sys.exit(1)

    i = int((x - blx) / dx)
    j = int((tly - y) / dy)
    print(x, y, i, j, field[j, i])

if __name__ == '__main__':
    main()
