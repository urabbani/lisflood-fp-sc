#!/usr/bin/env python3
import sys

lines = sys.stdin.readlines()

start = [i for (i, line) in enumerate(lines) if line.startswith("Stage information")]
end = [i for (i, line) in enumerate(lines) if line.startswith("Output, depths")]

start = start[0]
end = end[0]

stage_count = end - start - 2

z = []
# read in z values
for line in lines[start+1:end-1]:
    tokens = line.split()
    z.append(float(tokens[3]))

for line in lines[end+2:]:
    tokens = line.split()
    print(tokens[0], end=',')
    for i, token in enumerate(tokens[1:]):
        print(z[i] + float(token), end='')
        if i < len(tokens)-2:
            print(',', end='')
    print()
