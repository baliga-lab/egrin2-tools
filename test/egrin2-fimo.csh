#!/bin/csh

setenv PATH /tools/bin:${PATH}
setenv BATCHNUM `printf "%03d" $SGE_TASK_ID`
#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_mharris
#$ -t 1-3
#$ -M mharris@systemsbiology.org
#$ -cwd
#$ -pe serial 1
#$ -l mem_free=32G

python egrin2-fimo.py --cache_dir cache --seqs_file Escherichia_coli_K12_NC_000913.2 --output_dir "eco-out-$BATCHNUM/fimo_out" --progs_dir progs --overwrite true --input_dir eco-out-$BATCHNUM
