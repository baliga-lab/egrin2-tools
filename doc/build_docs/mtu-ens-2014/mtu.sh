#!/bin/csh

setenv LD_LIBRARY_PATH /tools/lib:/tools/R-3.0.3/lib64/R/lib
setenv PATH /tools/bin:${PATH}
setenv BATCHNUM `printf "%03d" $SGE_TASK_ID`
#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_abrooks
#$ -t 1-5
#$ -tc 8
#$ -l hostname="baliga2|baliga3"
#$ -M abrooks@systemsbiology.org
#$ -cwd
#$ -pe serial 20
#$ -l mem_free=10G

python cmonkey.py --organism mtu --ratios mtu-ens-2014/ratios-$BATCHNUM.tsv --config mtu-ens-2014/config-$BATCHNUM.ini --out mtu-out-$BATCHNUM --minimize_io

bzip2 -f mtu-out-$BATCHNUM/*.pkl
