#!/bin/bash

HELP="
Full pipeline to extract subtelomeres as FASTA file from assembly

Usage: extract [fasta-file, -x] [Options]

Options:
-n len        Extract subtelomeric sequence <len>-long
              default: 500000
-b buffer     Length of buffer for extraction
              should be longer than length of telomere
              setting shorter buffer will lead to errors
              default: 100000
-c #chrom     Only use first #chrom number of sequences from fasta-file
              default: 46
-t temp_dir   Directory for intermediate results
              default: temp
-o file       Where to put the output file
              the directory of the file must exist
              default: output.fa
-x            Assume the buffered subtelomeric sequence already extracted
              also needs the chromosome.sizes file
              this helps speedup subsequent reruns significantly
-d            delete temp files after run
              default:
-f bedfile    Also extract subtelomeres from bed file
              filename of the extracted feature file may be specified by option -B
-B f_file     Extract features into file f_file
              default: features.bed


Dependences:
faSize  conda   bioconda::ucsc-fasize
seqkit  conda   seqkit
seqtk   docker  staphb/seqtk
bigBedToBed docker  bioconda::ucsc-bigbedtobed
" # TODO

subtelo_len=500000
subtelo_buffer=100000
chrom_count=46
temp_dir="temp"
output_file="output.fa"
extract_subtelomere="true"
delete_temp=""
bed_file=""
subtelo_features="features.bed"

if [ $# -eq 0 ]; then
    echo "$HELP"
    exit
fi

if [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "$HELP"
    exit
fi

set -euo pipefail

fasta_file=$1
if [ ! "${fasta_file:0:1}" == "-" ]; then
    file=$1
    if [ "${file##*.}" == "gz" ]; then
        file="zcat $file"
    else
        file="cat $file"
    fi
    shift 1
fi

scriptdir=$(dirname $0)

while getopts "n:b:c:t:o:xf:B:d" o; do
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
        d)
            delete_temp="True"
            ;;
        f)
            bed_file=$OPTARG
            ;;
        B)
            subtelo_features=$OPTARG
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

    echo "Generating file for chromosome extraction"
    # Extract 500Kb + buffer as bed file
    ./$scriptdir/subtelomeres_bed.sh $temp_dir/chrom.sizes -o $temp_dir/seq_ends.bed -n $(( $subtelo_len + $subtelo_buffer )) -c $chrom_count

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

# TODO: reverse direction of subtelo_cropped
if [ $bed_file ]; then
    echo "Extracting from bedfile sequence"
    ./$scriptdir/bed_shifter.py \
        <(./$scriptdir/bed_shifter.py \
            <( bedtools intersect -a $bed_file -b $temp_dir/seq_ends.bed -wa ) \
            $temp_dir/seq_ends.bed \
        ) \
        <(./$scriptdir/rotate_end.py $temp_dir/subtelo_cropped.bed ) \
        > $subtelo_features
fi

cp $temp_dir/subtelo_cropped_t2c.fa $output_file

if [ -n "$delete_temp" ]; then
    rm -rf $temp_dir
fi
