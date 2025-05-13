#!/bin/bash

HELP="""
Adds splitter N-sequences and merges sequences
Also allows for filtering sequences

Usage: split-merge file.fa [Options]

Options:
-n split-len        Split using sequence of split-len N's
                    default 50000
-c char             Use char instead of N for splitting
-o out-file         Write output to out-file
                    default split-merge.fa
-t temp-dir         Temp dir
                    default temp-split-merge
-p pattern          Grep file.fa on name with this pattern
                    regex patterns are allowed
                    By default all sequences are added
-r rename           Pattern to use with sed
                    this is useful for correction in sorting
                    Script always adds 0's to >chr[num] at the beginning of name
                    Modification using sed are made after the grep phase
                    Use another separator such as pipe '|' instead of usual '/'
"""


if [ $# -eq 0 ]; then
    echo "$HELP"
    exit
fi

if [ $1 == "-h" ]; then
    echo "$HELP"
    exit
fi

files="$1"
shift 1

split=50000
split_char="N"
output="split-merge.fa"
temp="temp-split-merge"
grep_pattern=""
sed_pattern="-e s/^>chr([0-9][^0-9])/>chr00\1/ -e s/^>chr([0-9][0-9][^0-9])/>chr0\1/"

while getopts "n:c:o:t:p:r:h" o; do
    case $o in
        n)
            split=$OPTARG
            ;;
        c)
            split_char=$OPTARG
            ;;
        o)
            output=$OPTARG
            ;;
        t)
            temp=$OPTARG
            ;;
        p)
            grep_pattern="$grep_pattern -p $OPTARG"
            ;;
        r)
            sed_pattern="$sed_pattern -e $OPTARG"
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

if [ -n "$grep_pattern" ]; then
    seqkit grep --by-name $files -r $grep_pattern  > $temp/grepped.fa
else
    cat $files > $temp/grepped.fa
fi

grep ">" $temp/grepped.fa | \
    while read seq; do
        echo $seq
        printf "${split_char}%.0s" $(seq $split)
        echo ""
    done > $temp/splitter.fa

echo ">merged" > $output
seqkit concat $temp/grepped.fa $temp/splitter.fa | \
    sed -r $sed_pattern | \
    seqkit sort --by-name | \
    tee >( grep ">" > $output.seqs ) | \
    seqkit seq --seq >> "$output"
