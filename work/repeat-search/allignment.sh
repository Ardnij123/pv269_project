#!/bin/bash

HELP="""
Extracts alligned subsequences from .bedpe file and allings them using Clustal

Usage: allginment file.fa file.bedpe [Options]

Options:
-t temp-dir     Temp dir
                default temp-allignment
"""


if [ $# -eq 0 ]; then
    echo "$HELP"
    exit
fi

if [ $1 == "-h" ]; then
    echo "$HELP"
    exit
fi


fasta="$1"
bedpe="$2"

temp="temp-allignment"

shift 2

while getopts "t:h" o; do
    case $o in
        t)
            temp=$OPTARG
            ;;
        h)
            echo "$HELP"
            exit
            ;;
        *)
            echo "Unrecognized parameter $o. Printing help."
            echo "$HELP"
            exit
            ;;
    esac
done

mkdir -p $temp

# TODO: this is in bedpe format, should sub 1 from start
cat \
    <( tail -n +2 $bedpe | cut -f 1-3 | uniq ) \
    <( tail -n +2 $bedpe | cut -f 4-6 | sort -u | uniq ) \
    | sort -u | uniq > $temp/sequences.bed

seqkit subseq --bed $temp/sequences.bed $fasta > $temp/sequences.fa


