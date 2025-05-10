#!/bin/bash

HELP="
Creates self-maps for all fasta files specified as parameter

Usage: ./self_map.sh 'files' [-i id-perc] [-r resolution] [-o diretory]

Options:
-i id-perc      At least id-perc percentage of identity is needed to plot
                default 80
-r resolution   Resolution of final graph
                default 1000
-o directory    Output directory
                default self-maps
-p max-procs    Maximal number of concurent processes (used processors)
                default 4
-t temp_dir     Temporal directory for splitting fasta files
                default temp-split
-h              Print this help and exit

Needs ModDotPlot venv activated
"

id_perc=80
res=1000
out_dir="self-maps"
processors=4
temp_dir="temp-split"

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

while getopts "i:r:o:p:t:h" o; do
    case $o in
        i)
            id_perc=$OPTARG
            ;;
        r)
            res=$OPTARG
            ;;
        o)
            out_dir=$OPTARG
            ;;
        p)
            processors=$OPTARG
            ;;
        t)
            temp_dir=$OPTARG
            ;;

        h)
            echo -e "$HELP"
            exit
            ;;
        *)
            echo "Unrecognized parameter $o, priniting help"
            echo -e "$HELP"
            exit
            ;;
    esac
done


seqkit split --by-id $files --out-dir "$temp_dir"

ls $temp_dir/* | \
    xargs -i --max-procs=$processors bash -c \
    "echo starting {}; moddotplot static -f {} -id $id_perc -r $res -o $out_dir --deraster > /dev/null; echo {} finished"
