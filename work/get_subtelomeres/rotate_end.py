#!/bin/python3.11

import argparse


HELP="""
Rotates bed file features containing subsequence 'END'.
This script is for use with bed-shifter in this directory.

Usage: rotate bedfile
"""

parser = argparse.ArgumentParser(
	prog='bed-shifter',
	description=HELP,
)

parser.add_argument('bedfile', help='Bedfile to transform')

args = parser.parse_args()

features = args.bedfile

with open(features, 'r') as f:
    for line in f:
        feature = line.split()
        _, start, end = feature[0], int(feature[1]), int(feature[2])

        if "END" in line:
            start, end = end-1, start-1

        print("\t".join([feature[0]] + list(map(str, [start, end])) + feature[3:]))
