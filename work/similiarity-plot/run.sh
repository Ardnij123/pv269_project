#!/bin/bash

HELP="
Pipeline to create visualisations and plots for sequence allignment using moddotplot
This pipeline does not run you out of memory as it computes the allignments one by one

Usage: similiarity fasta-files [Options]

Options:
-i id-perc      At least id-perc percentage of identity is needed to plot
                default 80
-r resolution   Resolution of final graph
                default 1000
-o directory    Output directory
                default sim-maps
-p max-procs    Maximal number of concurent processes (used processors)
                default 4
-t temp_dir     Temporal directory for splitting fasta files
                default temp-split
-s value        Create self-plots
                default True
-c              Create cross-plots
                default False
-g              Create grid from outputted dotplots
                This option will also turn on cross-plots generation
                default False
-G file         Choose destination for outputted grid
                default dotplot-grid.png
-h              Print this help and exit

Needs ModDotPlot venv activated
"


scriptdir=$(dirname $0)

id_perc=80
res=1000
out_dir="sim-maps"
processors=4
temp_dir="temp-split"
self_plots="True"
cross_plots=""

grid_output=""
grid_file="dotplot-grid.png"

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

while getopts "i:r:o:p:t:hs:cgG:" o; do
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
        s)
            self_plots=$OPTARG
            ;;
        c)
            cross_plots="True"
            ;;
        h)
            echo -e "$HELP"
            exit
            ;;
        g)
            grid_output="True"
            cross_plots="True"
            ;;
        G)
            grid_file=$OPTARG
            ;;
        *)
            echo "Unrecognized parameter $o, priniting help"
            echo -e "$HELP"
            exit
            ;;
    esac
done


mdp_options="-i $id_perc -r $res -o $out_dir"

seqkit split --by-id $files --out-dir "$temp_dir" 2> /dev/null

if [ -n "$self_plots" ]; then
    echo "Plotting self-plots"

    ls $temp_dir/* | \
        xargs -i --max-procs=$processors bash -c \
        "moddotplot static -f {} $mdp_options --deraster > /dev/null"
fi

if [ -n "$cross_plots" ]; then
    echo "Plotting cross-plots"
    ls $temp_dir/* | \
        xargs -I'{1}' --max-procs=1 bash -c "echo 'Alligning file {1}';
            ls $temp_dir/* | xargs -I'{2}' --max-procs=$processors bash -c \
            \"moddotplot static -f {1} {2} $mdp_options --compare-only --width 1 --dpi 500 > /dev/null 2> /dev/null\";
            echo Done; rm $out_dir/*_COMPARE.svg; rm $out_dir/*_COMPARE_HIST.svg
        "
fi

if [ -n $grid_output ]; then
    ./$scriptdir/grid_plots.sh $files $out_dir
fi
