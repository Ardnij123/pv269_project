#!/bin/python3.11

import argparse


HELP="""
Shifts bed file features to correspond with another coordinate system.

Usage: bedshift bedfile shift

Example:
Assume a feature in bedfile:
chrA    235 584     ...

And a shift in shift-system (0-based, non-inclusive end, same as bed format):
chrA    48  2391    chrZ

This script will shift the feature 48 bases to the left as the new coordinate
system chrZ starts at base 48 of old coordinate system chrA.
Outputted file will have a single feature:
chrA    187 536     ...

There may be multiple coordinate systems deriving from a coordinate system.
In that case, the feature will be transposed into each of them.
If the transposition would result into being out of borders, the transposition
is thrown away. Should the intersection of feature with the new coordinate
system be only partial, the intersection is added as a new coordinates of feature.

Script may even transform the direction of the coordinate system. In that case
the coordinate line in shift-system should have higher number in begin
than in end. Again, the first number (this time higher) is inclusive while second
(lower) number is noninclusive. (To get mapping for same region, swithch
the numbers and add or subtract 1.)
Such mapping will switch the orientation of strand between '+' and '-'
"""

parser = argparse.ArgumentParser(
	prog='bed-shifter',
	description=HELP,
)

parser.add_argument('bedfile', help='Bedfile to transform')
parser.add_argument('shift', help='Shift system')

args = parser.parse_args()

features = args.bedfile
shifts = args.shift

shift_d = dict()

with open(shifts, 'r') as f:
    for line in f:
        system = line.split()
        if len(system) == 3:
            old, start, end = system
            new = old
        elif len(system) >= 4:
            old, start, end, new = system[:4]
        else:
            print("Error on parsing shift file: ", system)
            exit(1)

        if old not in shift_d:
            shift_d[old] = []
        shift_d[old].append((int(start), int(end), new))

with open(features, 'r') as f:
    for line in f:
        feature = line.split()
        old, ostart, oend = feature[0], int(feature[1]), int(feature[2])

        systems = shift_d.get(old, [])
        for start, end, new in systems:
            _ostart, _oend = ostart, oend
            reverse = (start >= end)

            if reverse:
                # new system is reversed -> reverse everything
                _ostart, _oend = _oend-1, _ostart-1
                start, end, _ostart, _oend = -start, -end, -_ostart, -_oend
            
            _ostart = max(start, _ostart)
            _oend = min(end, _oend)

            if _ostart >= _oend:
                continue

            _ostart -= start
            _oend -= start

            if reverse and len(feature) >= 6:
                strand = feature[5]
                if strand == '+':
                    strand = '-'
                elif strand == '-':
                    strand = '+'
                print("\t".join(
                    list(map(str, [new, _ostart, _oend]))
                    + feature[3:5] + [strand] + feature[6:]))
            else:
                print("\t".join(list(map(str, [new, _ostart, _oend])) + feature[3:]))
