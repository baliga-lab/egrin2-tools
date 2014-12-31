#!/usr/bin/env python
"""cMonkeyQSub.py
This is a tool to prepare a directory for ensemble runs. Given an organism code
and the ratios, matrix, split up the matrix into sub matrices and create
SGE qsub shell scripts

Example:

With pre-defined condition blocks:

python cMonkeyQSub.py --organism mtb --ratios 20141130.MTB.all.ratios.csv --targetdir mtb-ens-20141230 --numruns 20 --blocks 20141202.MTB.EGRIN2.blocks.csv --inclusion 20141202.MTB.EGRIN2.exclusion.blocks.csv --exclusion 20141202.MTB.EGRIN2.inclusion.blocks.csv --pipeline setenrich_pipeline.json --setenrich chipseq,tfoe --setenrich_files 20140725.MTB.ChIPSeq.csv,20140725.MTB.DE.csv

Without pre-defined condition blocks:

python cMonkeyQSub.py --organism mtb --ratios 20141130.MTB.all.ratios.csv --targetdir mtb-ens-20141230 --numruns 20 --mincols 50 --num_cores 1
"""
import argparse
import os
import cmonkey.datamatrix as dm
from egrin2.cMonkeyIniGen import *
from egrin2.ensemblePicker import * 

DESCRIPTION = """ensemble.py - prepare cluster runs"""

# Templates for Bourne Shell
QSUB_TEMPLATE_HEADER = """#!/bin/bash

export LD_LIBRARY_PATH=/tools/lib:/tools/R-3.0.3/lib64/R/lib
export PATH=/tools/bin:${PATH}
export BATCHNUM=`printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE = """#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -t 1-%d
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

python cmonkey.py --organism %s --ratios %s --config %s --out %s --minmize_io %s

bzip2 -f %s/*.pkl
"""

# Templates for csh

QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh

setenv LD_LIBRARY_PATH /tools/lib:/tools/R-3.0.3/lib64/R/lib
setenv PATH /tools/bin:${PATH}
setenv BATCHNUM `printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE_CSH = """#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -t 1-%d
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

python cmonkey.py --organism %s --ratios %s --config %s --out %s --minmize_io %s

bzip2 -f %s/*.pkl
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="Organism code")
    parser.add_argument('--ratios', required=True)
    parser.add_argument('--targetdir', required=True)
    parser.add_argument('--numruns', type=int, default=4)
    parser.add_argument('--mincols', type=int, default=8)
    parser.add_argument('--num_cores', type=int, default=1)
    parser.add_argument('--user', default=None)
    parser.add_argument('--csh', action='store_true')
    parser.add_argument('--blocks', default=None)
    parser.add_argument('--inclusion', default=None)
    parser.add_argument('--exclusion', default=None)
    parser.add_argument('--pipeline', default=None)
    parser.add_argument('--setenrich', default=None)
    parser.add_argument('--setenrich_files', default=None)
    args = parser.parse_args()

    if args.csh:
        header = QSUB_TEMPLATE_HEADER_CSH
        template = QSUB_TEMPLATE_CSH
    else:
        header = QSUB_TEMPLATE_HEADER
        template = QSUB_TEMPLATE

    # write ratios files
    if args.blocks is None:
      # if inclusion/exlcusion blocks are not defined, simply choose at random
      dm.prepare_ensemble_matrix( args.ratios, args.targetdir, args.numfiles,
                                 args.mincols )
    else:
      cols = ensemblePicker( blocks = args.blocks, inclusion = args.inclusion, exclusion = args.exclusion, nruns = args.numruns, ratios_file = args.targetdir )
      cols.pickCols_all()

    # write config files
    for i in range( 1, args.nruns+1 ):
      ini = cMonkeyIniGen( )
      ini.writeIni( os.path.join( args.targetdir, "config-%03d.ini" % i ) )

    with open(os.path.join(args.targetdir, "%s.sh" % args.organism), 'w') as outfile:
        if args.user is not None:
            login = args.user
        else:
            login = os.getlogin()

        if args.pipeline is not None:
          pipeline = "--pipeline " + args.pipeline
        else:
          pipeline = ""

        outfile.write(header)
        outfile.write(template % (login, args.numruns, login,
                                  args.num_cores,
                                  args.organism,
                                  os.path.join(args.targetdir, "ratios-$BATCHNUM.tsv"),
                                  os.path.join(args.targetdir, "config-$BATCHNUM.ini"),
                                  "%s-out-$BATCHNUM" % (args.organism),
                                  pipeline ) )
