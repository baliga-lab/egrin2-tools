#!/usr/bin/env python3
"""generate_cm2_runs.py
This is a tool to prepare a directory for ensemble runs. Given an organism code
and the ratios, matrix, split up the matrix into sub matrices and create
SGE qsub shell scripts

Example:

With pre-defined condition blocks:

python cm2sge.py --organism mtu --ratios 20141130.MTB.all.ratios.csv --targetdir mtu-ens-20140120 --numruns 500 --blocks 20141202.MTB.EGRIN2.blocks.csv --inclusion 20141202.MTB.EGRIN2.inclusion.blocks.csv --exclusion 20141202.MTB.EGRIN2.exclusion.blocks.csv --pipeline setenrich_pipeline.json --setenrich chipseq,tfoe --setenrich_files 20140725.MTB.ChIPSeq.csv,20140725.MTB.DE.csv --csh

Without pre-defined condition blocks:

python generate_cm2_runs.py --organism mtu --ratios 20141130.MTB.all.ratios.csv --targetdir mtb-ens-20141230 --numruns 20 --mincols 50 --num_cores 1 --csh
"""
import argparse
import os
import itertools
import logging
import random
import shutil

import ensemble.datamatrix as dm
import ensemble.cmconfig as cmconfig
import ensemble.ensemble as ensemble


DESCRIPTION = "generate_cm2_runs.py - prepare cluster runs for Sun Grid Engine"

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
#$ -tc %d
#$ -l hostname="baliga2|baliga3"
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=10G

cmonkey2 --organism %s --config %s --out %s --minimize_io %s %s

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
#$ -tc %d
#$ -l hostname="baliga2|baliga3"
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=10G

cmonkey2 --organism %s --config %s --out %s --minimize_io %s %s

bzip2 -f %s/*.pkl
"""

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None  # "make_ensemble.log"


def main():
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="3 letter organism code")
    parser.add_argument('--ratios', required=True, help="Path to ratios file")
    parser.add_argument('--targetdir', required=True, help="Path to output directory")
    parser.add_argument('--numruns', type=int, default=100, help="Number of cMonkey2 runs to configure")
    parser.add_argument('--ncbi_code', default="", help="NCBI organism code")
    parser.add_argument('--mincols', type=int, default=100, help="Minimum number of experiments to include in a cMonkey2 run")
    parser.add_argument('--num_cores', type=int, default=8,  help="Number of cores on cluster to request")
    parser.add_argument('--max_tasks', type=int, default=20, help="Maximum number of jobs to be sent to the cluster at a time")
    parser.add_argument('--user', default=None, help="Cluster user name")
    parser.add_argument('--csh', action='store_true', help="Flag to indicate C Shell")
    parser.add_argument('--blocks', default=None, help="Path to block definitions")
    parser.add_argument('--inclusion', default=None, help="Path to inclusion block definitions")
    parser.add_argument('--exclusion', default=None, help="Path to exclusion block definitions")
    parser.add_argument('--pipeline', default=None, help="Path to scoring pipeline config file")
    parser.add_argument('--setenrich', default=None, help="Name(s) of set enrichment 'sets' to include. Names should be comma separated.")
    parser.add_argument('--setenrich_files', default=None, help="Set enrichment files. File paths should be comma separated.")
    parser.add_argument('--rsat_base_url', default=None, help="Alternative RSAT base URL.")
    args = parser.parse_args()

    if args.csh:
        header = QSUB_TEMPLATE_HEADER_CSH
        template = QSUB_TEMPLATE_CSH
    else:
        header = QSUB_TEMPLATE_HEADER
        template = QSUB_TEMPLATE

    if not os.path.exists(args.targetdir):
        os.makedirs(args.targetdir)

    # write ratios files
    logging.info("Choosing ensemble conditions")
    if args.blocks is None:
        # if inclusion/exclusion blocks are not defined, simply choose at random
        logging.info("no inclusion/exclusion blocks defined, performing random sub matrix generation...")
        dm.prepare_ensemble_matrix(args.ratios, args.targetdir, args.numruns,
                                   args.mincols)
    else:
        ensemble.make_ensemble_ratios(args.ratios, args.blocks, args.exclusion, args.inclusion,
                                      args.numruns, args.targetdir)

    enrich_sets = args.setenrich.split(",")
    setenrich_files = args.setenrich_files.split(",")
    stripped_setenrich_files = [f.split('/')[-1] for f in setenrich_files]
    cmconfig.make_config_files(args.num_cores, enrich_sets, stripped_setenrich_files,
                               args.numruns, args.pipeline.split('/')[-1], args.targetdir)
    # copy auxiliary files
    shutil.copy(args.pipeline, args.targetdir)
    for f in setenrich_files:
        shutil.copy(f, args.targetdir)

    with open(os.path.join(args.targetdir, "%s.sh" % args.organism), 'w') as outfile:
        login = args.user if args.user is not None else os.getlogin()
        outfile.write(header)
        rsat_base_url = ""
        if args.rsat_base_url is not None:
            rsat_base_url = "--rsat_base_url %s" % args.rsat_base_url
        outfile.write(template % (login,
                                  args.numruns,
                                  args.num_cores,
                                  login,
                                  args.max_tasks,
                                  args.organism,
                                  "config-$BATCHNUM.ini",
                                  "%s-out-$BATCHNUM" % (args.organism),
                                  rsat_base_url,
                                  "ratios-$BATCHNUM.tsv",
                                  "%s-out-$BATCHNUM" % (args.organism)))
        logging.info("Done")
