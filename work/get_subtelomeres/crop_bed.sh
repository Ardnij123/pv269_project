#!/bin/bash

HELP="
Crops sequences given by BED file to some length

The cropping is given by the end of chromosome field of BED file.
By default '..._START' and '..._END' should be used

Usage: crop bed_file [-n len] [-s start] [-e end]

Options:
-n len     length of subsequence to leave
           default value '500000'
-s start   set different string as start
           default value '_START'
-e end     set different string as end
           default value '_END'
"

len=500000
start='_START'
end='_END'

if [ $# -eq 0 ]; then
    echo -e "$HELP"
    exit
fi

file=$1
shift 1

if [ ! -e $file ]; then
    echo "File $file does not exist."
    exit
fi

while getopts "n:s:e:" o; do
    case "$o" in
        n)
            len=${OPTARG}
            ;;
        s)
            start=${OPTARG}
            ;;
        e)
            end=${OPTARG}
            ;;
        *)
            echo "Unrecognized parameter $o, printing help"
            echo -e "$HELP"
            exit
            ;;
    esac
done

cat "$file" | while read line; do
    chrom="`echo "$line" | cut -f1`"
    chrom_start="`echo "$line" | cut -f2`"
    chrom_end="`echo "$line" | cut -f3`"
    if [ -n "`echo $chrom | egrep "$start$"`" ]; then
        echo -e "$chrom\t$chrom_start\t$(($chrom_start + $len))"
    else
        echo -e "$chrom\t$(($chrom_end - $len))\t$chrom_end"
    fi
done
