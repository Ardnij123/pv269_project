#!/bin/bash
#
# Creates a grid from directory of png files
#
# Usage: grid fasta-file output-directory

fasta_file=$1
shift 1

simetry_dir=$1
shift 1

temp_dir="temp-grid"
output_file="full_grid.png"

# TODO: This discriminates X/Y chromosomes as lowest value instead of highest
sequences=$(\
    grep ">" $fasta_file --no-filename |
    sed 's/^.//' |
    sed 's/chr//' | sed 's/^X/9999X/' | sed 's/^Y/9999Y/'|
    sort -n |
    sed 's/^9999//' | sed 's/^/chr/'
)

mkdir $temp_dir

echo "$sequences" | while read seq; do
    convert -append `echo "$sequences" | while read seq2; do echo "$simetry_dir/${seq}_${seq2}_COMPARE.png"; done` "$temp_dir/row_$seq.png"
done

convert +append `echo "$sequences" | while read seq; do echo "$temp_dir/row_$seq.png"; done | tac` "$output_file"

rm -r $temp_dir
