#!/bin/bash
#
# Download data needed

mkdir data

# Download HG002 genome
wget "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz" > "data/hg002v1.1.fasta.gz"

# Download annotation data for methylation
wget "https://public.gi.ucsc.edu/~mcechova/HG002/Q100_HiFi_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed" > "data/Q100_HiFi_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed"
wget "https://public.gi.ucsc.edu/~mcechova/HG002/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5hmC.bed" > "data/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5hmC.bed"
wget "https://public.gi.ucsc.edu/~mcechova/HG002/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed" > "data/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed"
