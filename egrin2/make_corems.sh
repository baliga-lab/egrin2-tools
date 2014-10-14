#!/bin/bash

# This shell script controls the entire corem computation process.
# It delegates cmonkey data processing to python scripts, and the respective
# backbone and other computation to R and C++ scripts

readonly NUM_CORES=4
readonly EDGE_FILE="backbone.edges"


# Run compute_tanimoto on n separate cores, join the individual results
# into a single tanimoto.out file
# Parameters:
# 1: number of genes
# 2: number of cores
# 3: path to edge file
# returns: nothing
run_tanimoto() {
    local num_genes=$1
    local num_cores=$2
    local edge_file=$3
    local lower
    local upper
    local i

    local step=`expr $num_genes / $num_cores`
    local interv=(`seq 0 $step $num_genes`)

    local num_elems=${#interv[@]}
    local last=`expr $num_elems - 1`

    if [ ${interv[$last]} -ne $num_genes ]; then
        interv[$last]=$num_genes
    fi

    local until=`expr $num_elems - 2`

    for i in `seq 0 $until`
    do
        lower=${interv[$i]}
        upper=${interv[`expr $i + 1`]}
        echo "compute_tanimoto $edge_file 0 $lower $upper"
        compute_tanimoto $edge_file 0 $lower $upper &
    done
    wait
    echo "Tanimoto processes finished, joining results..."
    echo "cat $edge_file.tanimoto_* > tanimoto.out"
    cat $edge_file.tanimoto_* > tanimoto.out    
    echo "rm $edge_file.tanimoto_*"
    rm $edge_file.tanimoto_*
}

main() {
  Rscript extract_backbone.Rscript --matrix $1
  num_lines=`wc -l $1 | sed -e "s/ $1//"`
  num_genes=`expr $num_lines - 1`

  # this will create the *.numid2name and *.wpairs files
  adjmat2wpairs $EDGE_FILE 0 0
  run_tanimoto $num_genes $NUM_CORES $EDGE_FILE
  echo "DONE !"
}

usage() {
  echo "usage: make_corems.sh <matrix-file>"
}
if [ "$#" -ne 1 ]; then
    usage
else
    main $1
fi
