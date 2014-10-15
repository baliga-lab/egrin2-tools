#!/bin/bash
set -e

# This shell script controls the entire corem computation process.
# It delegates cmonkey data processing to python scripts, and the respective
# backbone and other computation to R and C++ scripts

readonly NUM_CORES=4
readonly EDGE_FILE="backbone.edges"
readonly INCREMENT=0.1

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
    echo "cat $edge_file.tanimoto_* > $edge_file.tanimoto"
    cat $edge_file.tanimoto_* > $edge_file.tanimoto
    echo "rm $edge_file.tanimoto_*"
    rm $edge_file.tanimoto_*
}

make_cluster_communities() {
    local edge_file=$1
    local i

    vals=`seq 0 $INCREMENT 1`
    for i in $vals
    do
        if [ "$i" != "0.0" ]; then
            echo "cluster_communities $edge_file $i"
            cluster_communities $edge_file $i &
        fi
    done

    # collect the results
    wait

    cat $edge_file.density_* > $edge_file.density
    echo "cat $edge_file.density_* > $edge_file.density"

    rm $edge_file.density_*
    echo "rm $edge_file.density_*"
}


main() {
  Rscript extract_backbone.Rscript --matrix $1
  num_lines=`wc -l $1 | sed -e "s/ $1//"`
  num_genes=`expr $num_lines - 1`

  # this will create the *.numid2name and *.wpairs files
  adjmat2wpairs $EDGE_FILE 0 0
  run_tanimoto $num_genes $NUM_CORES $EDGE_FILE
  make_cluster_communities $EDGE_FILE
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
