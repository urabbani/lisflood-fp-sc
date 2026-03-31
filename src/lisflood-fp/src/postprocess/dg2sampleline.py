#!/usr/bin/env python3
from osgeo import gdal
import math
import sys

def usage():
    print('Usage: dg2sampleline.py <sampleline_const> <sampleline_slope>')
    sys.exit(1)

def delta(const_lines):
    x1 = float(const_lines[0].split(' ')[2])
    y1 = float(const_lines[0].split(' ')[3])
    x2 = float(const_lines[1].split(' ')[2])
    y2 = float(const_lines[1].split(' ')[3])

    dx = x2 - x1
    dy = y2 - y1

    return dx, dy

def main():
    if len(sys.argv) <= 2:
        usage()

    with open(sys.argv[1]) as const, open(sys.argv[2]) as slope:
        const_lines = const.readlines()
        slope_lines = slope.readlines()

        assert len(const_lines) == len(slope_lines)

        dx, dy = delta(const_lines)

        for const_line, slope_line in zip(const_lines, slope_lines):
            [_, _, x_centre, y_centre, f0] = const_line.split(' ')
            x_centre = float(x_centre)
            y_centre = float(y_centre)
            f0 = float(f0)
            f1 = float(slope_line.split(' ')[4])

            limit_pos = f0 - math.sqrt(3.0)*f1
            limit_neg = f0 + math.sqrt(3.0)*f1

            print(x_centre-0.5*dx, y_centre-0.5*dy, limit_pos)
            print(x_centre+0.5*dx, y_centre+0.5*dy, limit_neg)

if __name__ == '__main__':
    main()

