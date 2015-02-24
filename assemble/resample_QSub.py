#!/usr/bin/env python
"""resample_QSub.py
This is the template for generating resample QSub scripts
for submission to SGE cluster

"""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import argparse
import os
import itertools

DESCRIPTION = """resample_QSub.py - generate QSub scripts for SGE cluster"""

# Templates for Bourne Shell
QSUB_TEMPLATE_HEADER = """#!/bin/bash

export PATH=/tools/bin:${PATH}
export BATCHNUM=`printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE = """#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -tc %d
#$ -l hostname="baliga1baliga2|baliga3"
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=10G

python cmonkey.py --organism %s --ratios %s --config %s --out %s

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
#$ -l hostname="baliga1baliga2|baliga3"
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=10G

python cmonkey.py --organism %s --ratios %s --config %s --out %s --minimize_io

bzip2 -f %s/*.pkl
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="Organism code")
    parser.add_argument('--ratios', required=True)
    parser.add_argument('--targetdir', required=True)
    parser.add_argument('--numruns', type=int, default=4)
    parser.add_argument('--ncbi_code', default="")
    parser.add_argument('--mincols', type=int, default=8)
    parser.add_argument('--num_cores', type=int, default=8)
    parser.add_argument('--max_tasks', type=int, default=20)
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


    with open(os.path.join(args.targetdir, "%s.sh" % args.organism), 'w') as outfile:
        if args.user is not None:
            login = args.user
        else:
            login = os.getlogin()

        outfile.write(header)
        outfile.write(template % (login, args.numruns, login,
                                  args.num_cores,
                                  args.max_tasks,
                                  args.organism,
                                  os.path.join(args.targetdir, "ratios-$BATCHNUM.tsv"),
                                  os.path.join(args.targetdir, "config-$BATCHNUM.ini"),
                                  "%s-out-$BATCHNUM" % (args.organism),
                                  "%s-out-$BATCHNUM" % (args.organism) ) )
