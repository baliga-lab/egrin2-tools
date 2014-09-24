#!/tools/bin/python

""" egrin2-coding_fracs.py
Author:  Micheleen M. Harris
Description:  Coding fractions step of the EGRIN2 pipeline.  This program makes a script which can be submitted to an SGE scheduler on a cluster.  It will calculate the fraction of motifs from FIMO that are in coding regions.

This script assumes an ensemble run under a directory structure that 
is as follows:

<dir>
  |-- <prefix>001
             |-- cmonkey_run.db
             |-- ...
  |-- <prefix>002
  |-- ...
  |-- <this is where a shell script will be run with qsub>

Instructions:  Run resulting script (<name you choose>.sh) under <dir> as shown above.  A script must be made for each cmonkey run.
Note:  if you are using c-shell make sure to use the --csh flag when running this program.

Example:  
python egrin2-coding_fracs.py --cache_dir cache --features_file cache/Escherichia_coli_K12_features --num_runs 3 --num_cores 1 --organism_name eco --user mharris --csh

Output will be in each <prefix> directory named coding_fracs.tsv.
"""

import optparse
import sys
import os
import glob
import csv

import time

# Templates for Bourne Shell
QSUB_TEMPLATE_HEADER = """#!/bin/bash

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
#$ -pe serial %s
#$ -l mem_free=32G

python run_coding_fracs.py --cache_dir %s --features_file %s --organism_name %s --input_dir %s"""


# Templates for csh

QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh

setenv PATH /tools/bin:${PATH}
set BATCHNUM="`printf '%03d' ${SGE_TASK_ID}`"
"""

QSUB_TEMPLATE_CSH = """#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -t 1-%d
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %s
#$ -l mem_free=32G

python run_coding_fracs.py --cache_dir %s --features_file %s --organism_name %s --input_dir %s"""

def main():

	#  Collect & check args

	op = optparse.OptionParser()
	op.add_option('-c', '--cache_dir', help='The cmpython cache directory where the genome features and info live')
	op.add_option('-f', '--features_file', help='The sequence file (with coding regions) (found in cache/<organism name>features')
	op.add_option('-m', '--num_runs', help='The number of cmonkey runs')
	op.add_option('-n', '--num_cores', help='Number of cores to use on cluster', default=1)
	op.add_option('-o', '--organism_name', help='The organism name (e.g. eco or hal)')
	op.add_option('-u', '--user', help='User name on cluster')
	op.add_option('-s', '--csh', help='If c-shell indicate with this flag', action='store_true')
	opt, args = op.parse_args()

	if not opt.cache_dir:
		op.error('need --cache_dir option.  Use -h for help.')
	if not opt.features_file:
		op.error('need --features_file option.  Use -h for help.')
	if not opt.organism_name:
		op.error('need --organism_name option.  Use -h for help.')
	if not opt.num_runs:
		op.error('need --num_runs option.  Use -h for help.')

	if opt.csh:
		header = QSUB_TEMPLATE_HEADER_CSH
		template = QSUB_TEMPLATE_CSH
	else:
		header = QSUB_TEMPLATE_HEADER
		template = QSUB_TEMPLATE

	with open(os.path.join(os.getcwd(), 'egrin2-coding_fracs.sh'), 'w') as outfile:
		if opt.user is not None:
			login = opt.user
		else:
			os.getlogin()

		outfile.write(header)
		outfile.write(template % (login, 
			int(opt.num_runs), 
			login,
			opt.num_cores, 
			opt.cache_dir, 
			opt.features_file,
			opt.organism_name,
			"%s-out-${BATCHNUM}" % (opt.organism_name)))


if __name__ == '__main__':
	main()


