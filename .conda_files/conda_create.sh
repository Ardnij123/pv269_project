#!/bin/bash

conda init
conda create -n pv269_project
conda install -n pv269_project seqkit
conda install -n pv269_project bioconda::ucsc-fasize
conda install -n pv269_project bioconda::mashmap
conda install -n pv269_project conda-forge::gnuplot
conda install -n pv269_project bioconda::bedtools
conda install -n pv269_project ucsc-bedtobigbed
conda install -n pv269_project ucsc-bigbedtobed
