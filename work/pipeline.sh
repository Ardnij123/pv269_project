#!/bin/bash
#
# Full pipeline for subtelomere extraction and analysis


scriptdir=$(dirname $0)

# TODO: init conda, moddotplot

# Extract subtelomeres
./$scriptdir/get_subtelomeres/run.sh $1 -d

# TODO: add mapping, visualisation, ...
