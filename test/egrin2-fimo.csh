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

python egrin2-fimo.py -g Seqs_file -o "eco-out-$BATCHNUM/fimo_out" -p progs -w true -i eco-out-$BATCHNUM
