#!/bin/bash
#
# Generates a bed file from fasta file with region spanning whole sequences
# of the fasta file

file=$1

faSize -detailed $file > ${file%.*}.sizes

cat ${file%.*}.sizes | while read seq; do
    id="`echo "$seq" | cut -f1`"
    length="`echo "$seq" | cut -f2`"
    echo -e "$id\t0\t$length"
done
