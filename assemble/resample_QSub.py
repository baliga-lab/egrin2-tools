#!/usr/bin/env python
"""resample_QSub.py
This is the template for generating resample QSub scripts
for submission to SGE cluster

python resample_QSub.py --name resample_55 --host primordial --db mtu_db --n_rows 55 --user abrooks

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

# Templates for csh

QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh

"""

QSUB_TEMPLATE_CSH = """#$ -S /bin/csh
#$ -N %(name)s
#$ -o 'out_messages.txt'
#$ -e 'error_messages.txt'
#$ -m be
#$ -q baliga
#$ -P Bal_%(user)s
#$ -l hostname="baliga2|baliga3"
#$ -M %(user)s@systemsbiology.org
#$ -cwd
#$ -pe serial %(n_cores)i
#$ -l mem_free=8G

python resample.py --host %(host)s --db %(db)s --n_rows %(n_rows)i --n_resamples %(n_resamples)i --port %(port)i

"""

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument( '--host', required=True, type=str, help="Host for MongoDB" )
    parser.add_argument('--db', required=True, type=str, help="Database name")
    parser.add_argument('--n_rows', required=True, type=int, help="Gene set size to test")
    parser.add_argument('--n_resamples', default=1000, type=int, help="Number of resamples to compute")
    parser.add_argument('--port', default=27017, help="MongoDB port", type=int )
    parser.add_argument('--cols', default=None, help="Columns (experiments) for resampling. Should be path to tab-delimited file containing column names that map to egrin2_col_names in MongoDB database, eg experiment names.", type=str )
    parser.add_argument('--user', default=None)
    parser.add_argument('--targetdir', default="./", help="Directory to store Qsub scripts")
    parser.add_argument('--n_cores', default=1, help="Number of cores to reserve on cluster")

    args = parser.parse_args()

    header = QSUB_TEMPLATE_HEADER_CSH
    template = QSUB_TEMPLATE_CSH

    args_d = vars(args)
    args_d["name"] = "r.%s" % args_d["n_rows"]

    with open(os.path.join(args.targetdir, "%s.sh" % args.name), 'w') as outfile:
        if args.user is not None:
            login = args.user
        else:
            login = os.getlogin()

        outfile.write(header)
        outfile.write(template % vars(args) )
