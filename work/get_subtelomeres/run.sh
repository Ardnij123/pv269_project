#!/bin/bash

HELP="
Full pipeline to extract subtelomeres as FASTA file from assembly

Usage: extract fasta-file [-n len] [-b buffer] [-c #chrom] [-t temp_dir] [-o file] [-x]

Options:
-n len        extract subtelomeric sequence <len>-long
              default: 500000
-b buffer     length of buffer for extraction
              should be longer than length of telomere
              setting shorter buffer will lead to errors
              default: 100000
-c #chrom     only use first #chrom number of sequences from fasta-file
              default: 46
-t temp_dir   directory for intermediate results
              default: temp
-o file       where to put the output file
              the directory of the file must exist
              default: output.fa
-x            assume the buffered subtelomeric sequence already extracted
              also needs the chromosome.sizes file
              this helps speedup subsequent reruns significantly


Dependences:
faSize  conda   bioconda::ucsc-fasize
seqkit  conda   seqkit
seqtk   docker  staphb/seqtk
bigBedToBed docker  bioconda::ucsc-bigbedtobed
"

subtelo_len=500000
subtelo_buffer=100000
chrom_count=46
temp_dir="temp"
output_file="output.fa"
extract_subtelomere="true"

if [ $# -eq 0 ]; then
    echo "$HELP"
    exit
fi

set -euo pipefail

fasta_file=$1
file=$1
if [ "${file##*.}" == "gz" ]; then
    file="zcat $file"
else
    file="cat $file"
fi
shift 1

scriptdir=$(dirname $0)

while getopts "n:b:c:t:o:x" o; do
    case $o in
        n)
            subtelo_len=$OPTARG
            ;;
        b)
            subtelo_buffer=$OPTARG
            ;;
        c)
            chrom_cout=$OPTARG
            ;;
        t)
            temp_dir=$OPTARG
            ;;
        o)
            output_file=$OPTARG
            ;;
        x)
            extract_subtelomere=""
            ;;
        *)
            echo "Unrecognized parameter $o, priniting help"
            echo -e "$HELP"
            exit
            ;;
    esac
done


mkdir -p $temp_dir

if [ $extract_subtelomere ]; then
    echo "Extracting chromosome sizes"
    # Extract chromosome sizes
    faSize -detailed <($file) > $temp_dir/chrom.sizes
fi

echo "Generating file for chromosome extraction"
# Extract 500Kb + buffer as bed file
./$scriptdir/subtelomeres_bed.sh $temp_dir/chrom.sizes -o $temp_dir/seq_ends.bed -n $(( $subtelo_len + $subtelo_buffer )) -c $chrom_count

if [ $extract_subtelomere ]; then
    echo "Extracting sequences close to ends of chromosomes"
    # Seqkit subseq for sequence ends
    # rename to >chr1_MATERNAL_START, >chr12_PATERNAL_END, ...
    seqkit subseq --quiet --bed $temp_dir/seq_ends.bed $fasta_file | \
    seqkit replace -p "(.)* " -r "" > $temp_dir/seq_ends.fa
fi

echo "Filtering out telomeric sequences"
# Filter out only telomeric subsequences
# Select telomeric sequences
seqtk="/home/jindmen/temp/seqtk/seqtk/seqtk"
$seqtk telo $temp_dir/seq_ends.fa > $temp_dir/telo.bed 2> /dev/null

# Invert to get subtelomeric sequences
# use bedtools to subtract $temp_dir/telo.bed from full sequences
./$scriptdir/full_bed.sh $temp_dir/seq_ends.fa > $temp_dir/seq_ends.full.bed
bedtools subtract -a $temp_dir/seq_ends.full.bed -b $temp_dir/telo.bed > $temp_dir/subtelo.bed

echo "Cropping subtelomeric sequences to exact size"
# Choose subsequence to only span 500Kb
./$scriptdir/crop_bed.sh $temp_dir/subtelo.bed -n $subtelo_len > $temp_dir/subtelo_cropped.bed

seqkit subseq --quiet --bed $temp_dir/subtelo_cropped.bed $temp_dir/seq_ends.fa |\
sed "s/_[^_]*$//" > $temp_dir/subtelo_cropped.fa

echo "Reversing end sequences"
seqkit concat -f \
    <(seqkit grep -nrp ".*_START" $temp_dir/subtelo_cropped.fa) \
    <(seqkit grep -nrp ".*_END" $temp_dir/subtelo_cropped.fa | \
        seqkit seq -rp -t"dna" -v) > $temp_dir/subtelo_cropped_t2c.fa 2> /dev/null

cp $temp_dir/subtelo_cropped_t2c.fa $output_file
