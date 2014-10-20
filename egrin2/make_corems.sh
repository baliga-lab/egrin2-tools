#!/bin/bash
set -e

# This shell script controls the entire corem computation process.
# It delegates cmonkey data processing to python scripts, and the respective
# backbone and other computation to R and C++ scripts

readonly NUM_PARTS=4
readonly EDGE_FILE="backbone.edges"
readonly INCREMENT=0.1

# Run compute_tanimoto on n separate cores, join the individual results
# into a single tanimoto.out file
# Parameters:
# 1: number of genes
# 2: number of partions to run tanimoto
# 3: path to edge file
# returns: nothing
run_tanimoto() {
    local num_genes=$1
    local num_parts=$2
    local edge_file=$3
    local i

    local step=`expr $num_genes / $num_parts`
    local interv=(`seq 0 $step $num_genes`)

    local num_elems=${#interv[@]}
    local last=`expr $num_elems - 1`

    # ensure all genes are included
    if [ ${interv[$last]} -ne $num_genes ]; then
        interv[$last]=$num_genes
    fi
    local pids=""
    for i in $(seq 0 $(expr $num_elems - 2))
    do
        echo "compute_tanimoto $edge_file 0 ${interv[$i]} ${interv[`expr $i + 1`]}"
        compute_tanimoto $edge_file 0 ${interv[$i]} ${interv[`expr $i + 1`]} &
        pids="$pids $!"
    done

    local num_errors=0
    for pid in $pids
    do
       echo "wait for job $pid..."
       wait $pid || let "num_errors+=1"
    done
    # we check for error
    if [ $num_errors -eq 0 ]; then
        echo "Tanimoto processes finished, joining results..."
        echo "cat $edge_file.tanimoto_* > $edge_file.tanimoto"
        cat $edge_file.tanimoto_* > $edge_file.tanimoto
        echo "rm $edge_file.tanimoto_*"
        rm $edge_file.tanimoto_*
        return 0
    else
        echo "errors while running compute_tanimoto"
        return 1
    fi
}

make_cluster_communities() {
    local edge_file=$1
    local i
    local num_errors=0
    local vals=`seq 0 $INCREMENT 1`
    local pids=""

    for i in $vals
    do
        if [ "$i" != "0.0" ]; then
            echo "cluster_communities $edge_file $i"
            cluster_communities $edge_file $i &
            pids="$pids $!"
        fi
    done

    for pid in $pids
    do
       echo "wait for job $pid..."
       wait $pid || let "num_errors+=1"
    done

    # collect the results
    if [ $num_errors -eq 0 ]; then
        cat $edge_file.density_* > $edge_file.density
        echo "cat $edge_file.density_* > $edge_file.density"

        rm $edge_file.density_*
        echo "rm $edge_file.density_*"
        return 0
    else
        echo "error while running cluster_communities"
        return 1
    fi
}

main() {
  Rscript extract_backbone.Rscript --matrix $1
  local num_genes=$(expr $(cat $1 | wc -l) - 1)
  local retcode

  # this will create the *.numid2name and *.wpairs files
  adjmat2wpairs $EDGE_FILE 0 0

  if ! run_tanimoto $num_genes $NUM_PARTS $EDGE_FILE ; then
      echo "ERROR: running compute_tanimoto incomplete"
      exit 1
  fi

  if ! make_cluster_communities $EDGE_FILE ; then
      echo "ERROR: running cluster_communities incomplete"
      exit 1
  fi

  cutoff=`Rscript choose_cutoff.Rscript --densityfile $EDGE_FILE.density`
  echo "DONE with cutoff: $cutoff"
  # the result of get_communities is a tab separated file with the 5 columns
  # gene1 gene2 community_id community_density community_weighted_density
  communities_file=$EDGE_FILE.communities_$cutoff
  echo "getting_communities -> $communities_file"
  getting_communities $EDGE_FILE $cutoff
  # filter corems where Community.Weighted.Density > 0
  cat $communities_file | ./filter_communities.py > $communities_file.filtered

  # TODO: load into filtered communities into gene-by-gene matrix
  # TODO: normalize ratios ?
}

usage() {
  echo "usage: make_corems.sh <matrix-file>"
}
if [ "$#" -ne 1 ]; then
    usage
else
    main $1
fi
