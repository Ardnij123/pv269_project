#!/bin/bash

HELP="
This script extracts subtelomere region sequences as a BED file from
a chromosome.sizes file. The script only takes chrom_number first/largest
sequneces from the fasta file into account.

Usage: extract sizes_file [-n len] [-c count] [-s] [-o output]
       extract -h

Options:
-n len      length of subtelomeres to take in base pairs
            default value '500000'
-c count    only take subtelomeres from first count sequences
            default value '24' (chr1-22, X, Y)
-s          sort the sizes file by length of subsequences
-h          print this help
-o output   choose different file for the output
            default value 'subtelomeres'
"

len=500000
count=24
sorting=0
output="subtelomeres"

if [ $# -eq 0 ]; then
    echo -e "$HELP"
    exit
fi

if [ "$1" == "-h" ]; then
    echo -e "$HELP"
    exit
fi

if [ ! -e "$1" ]; then
    echo "File $1 does not exist."
    exit
fi

file=$1
shift 1

while getopts "n:c:so:h" o; do
	case "${o}" in
		n)
			len=${OPTARG}
			;;
		c)
			count=${OPTARG}
			;;
		s)
            sorting=1
			;;
        o)
            output=${OPTARG}
            ;;
        h)
            echo -e "$HELP"
            exit
            ;;
        *)
            echo "Unrecognized parameter ${o}, printing help"
            echo -e "$HELP"
            exit
            ;;
	esac
done

head -n "$count" "$file" | while read line; do
    chrom="`echo "$line" | cut -f1`"
    size="`echo "$line" | cut -f2`"
    echo -e "$chrom\t0\t$len\t${chrom}_START"
    echo -e "$chrom\t$(( $size - $len ))\t$size\t${chrom}_END"
done > "$output"
