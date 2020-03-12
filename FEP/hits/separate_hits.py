#!/usr/bin/env python

import sys

with open('hits.sdf') as f:
    lines = f.readlines()
lines = [x.split('\n')[0] for x in lines]
hit = 1
for x,line in enumerate(lines):
    if 'RDKit' in line:
        start = x
    if '$$$$' in line:
        if x != start:
            with open('hit%d.sdf'%hit,'w') as file:
                for n in range(start,x+1):
                    file.write(lines[n])
        hit += 1
