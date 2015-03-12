#!/tools/bin/python

import optparse
import sys
import os
import glob
import csv

import time
import bz2

import pandas as pd
import numpy as np
import numpy.core.defchararray as npstr

""" coding_fracs.py
Author:  Micheleen M. Harris, David J Reiss
Description:  Coding fractions step of the EGRIN2 pipeline.  This program makes a script which can be 
submitted to an SGE scheduler on a cluster.  It will calculate the fraction of motifs 
that are in coding regions based on FIMO scans.

This script assumes an ensemble run under a directory structure that 
is as follows:

<dir>
  |-- <prefix>001
             |-- cmonkey_run.db
             |-- fimo-outs/fimo-out*
             |-- ...
  |-- <prefix>002
  |-- ...
  |-- <this is where a shell script will be run with qsub>

and that the fimo.py script has been run and the fimo_batch_script.sh has been submitted to the
cluster (and completed).

Instructions:  Run resulting script (<name you choose>.sh) under <dir> as shown above.  
A script must be made for each cmonkey run.
Note:  if you are using c-shell make sure to use the --csh flag when running this program.

Example:  
python coding_fracs.py --features cache/Escherichia_coli_K12_features
                       --organism eco --user mharris --csh

Output will be in each <prefix> directory named coding_fracs.tsv.
"""

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
#$ -l mem_free=8G

python egrin2-tools/postproc/coding_fracs.py --features %s --organism %s --input_dir %s > %s/coding_fracs.out
"""


# Templates for csh

QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh -f

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
#$ -l mem_free=8G

python egrin2-tools/postproc/coding_fracs.py --features %s --organism %s --input_dir %s >& %s/coding_fracs.out
"""

def get_in_coding_rgn( input_dir, features ):
    fimo_files = np.sort( np.array( glob.glob(os.path.join(input_dir, "fimo-outs/fimo-out-*")) ) )
    print 'Number of fimo files = ', str(len(fimo_files))

    #  Get coding regions from features file (e.g. Escherichia_coli_K12_features)
    f = open(features, 'r')
    skip = 1
    line = f.readline()
    while 'header' not in line:
        line = f.readline()
        skip += 1
    f.close()

    features = pd.read_table(features,skiprows=skip) ## engine='python',
    features = features[ features.type != 'SEQ_END' ]
    start_pos = npstr.replace(features.start_pos.values.astype(str),'<','').astype(int)
    end_pos = npstr.replace(features.end_pos.values.astype(str),'>','').astype(int)

    total_coding_fracs = {}
    for f in fimo_files:
        print f
        fimo = pd.read_table( bz2.BZ2File(f), sep='\t' )
        if np.in1d('in_coding_rgn', fimo.columns)[0]:
            print 'SKIPPING', f, '; already done'
            #total_coding_fracs.append( np.nanmean( fimo.in_coding_rgn ) )
            #continue
        else:
            is_bad = np.zeros( fimo.shape[0], dtype=bool )
            for i in xrange(fimo.shape[0]):
                row = fimo.ix[i]
                hits = np.sum( (start_pos <= row.start) & (end_pos >= row.start) )
                if hits <= 0:
                    hits = np.sum( (start_pos <= row.stop) & (end_pos >= row.stop) )
                if hits > 0:
                    is_bad[i] = True
            
            ## write out fimo file with new column
            fimo['in_coding_rgn'] = is_bad
            fimo.to_csv( bz2.BZ2File( f, 'w' ), sep='\t', index=False )

        if fimo.shape[0] <= 0:
            continue
        ## write out summary for each motif in the run
        grpd = fimo.groupby('#pattern name').mean()
        mean_is_bad = grpd['in_coding_rgn'].values
        mot_ind = grpd.index.values
        mean_is_bad[ mean_is_bad == True ] = 1.0
        mean_is_bad[ mean_is_bad == False ] = 0.0
        ff = os.path.basename(f).replace('fimo-out-','').replace('.bz2','')
        for i in range(len(mean_is_bad)):
            total_coding_fracs[ff+'_'+str(mot_ind[i])] = round(mean_is_bad[i], 4)
    
    coding_fracs = {'motif':[a for a in sorted(total_coding_fracs.keys())], \
                        'coding_frac':[total_coding_fracs[a] for a in sorted(total_coding_fracs.keys())]}
    coding_fracs = pd.DataFrame(coding_fracs)
    coding_fracs.to_csv( bz2.BZ2File( os.path.join(input_dir, "coding_fracs.tsv.bz2" ), 'w' ), 
                         sep='\t', index=False )
    return coding_fracs


def get_total_coding_rgn( features ):
    """get total fraction of genome that is covered by coding region. 
       Not currently used in EGRIN2 pipeline but can be used to assess validity of cis-regulatory motifs."""
    #  Get coding regions from features file (e.g. Escherichia_coli_K12_features)
    f = open(features, 'r')
    skip = 1
    line = f.readline()
    while 'header' not in line:
        line = f.readline()
        skip += 1
    f.close()

    features = pd.read_table(features,skiprows=skip)
    genome_len = int(features.ix[0].end_pos)
    features = features[ features.type != 'SEQ_END' ]

    start_pos = npstr.replace(features.start_pos.values.astype(str),'<','').astype(int)
    end_pos = npstr.replace(features.end_pos.values.astype(str),'>','').astype(int)

    hits = np.zeros( genome_len + 1, dtype=bool ) ## tbd: use bitarray? (package bitarray https://pypi.python.org/pypi/bitarray)
    for i in xrange(len(start_pos)):
        hits[ start_pos[i]:end_pos[i] ] = True

    return np.mean( hits )

def main():
    #  Collect & check args

    op = optparse.OptionParser()
    ##op.add_option('-c', '--cache_dir', default='cache', help='The cmpython cache directory where the genome features and info live')
    op.add_option('-f', '--features', help='The features file (with coding regions) (found in cache/<organism name>features')
    ##op.add_option('-m', '--num_runs', default=100, help='The number of cmonkey runs')
    ##op.add_option('-n', '--num_cores', default=2, help='Number of cores to use on cluster')
    op.add_option('-o', '--organism', help='The organism name (e.g. eco or hal)')
    op.add_option('-u', '--user', help='User name on cluster')
    op.add_option('-s', '--csh', help='If c-shell indicate with this flag', action='store_true')

    ## option for actually doing the computation. if this is None then we're only writing out the qsub script
    op.add_option('-i', '--input_dir', help='The cmonkey run input dir where the fimo files live (e.g. eco-out-001)')

    opt, args = op.parse_args()

    #if not opt.cache_dir:
    #    op.error('need --cache_dir option.  Use -h for help.')
    if not opt.features:
        op.error('need --features option.  Use -h for help.')
    if not opt.organism:
        op.error('need --organism option.  Use -h for help.')
    ##if not opt.num_runs:
    ##    op.error('need --num_runs option.  Use -h for help.')

    org_out_dirs = glob.glob("%s-out-*" % opt.organism)

    if opt.input_dir:
        ## do the coding rgns calc, append column to fimo files
        get_in_coding_rgn( opt.input_dir, opt.features )

    else:
        ## write out qsub script
        if opt.csh:
            header = QSUB_TEMPLATE_HEADER_CSH
            template = QSUB_TEMPLATE_CSH
        else:
            header = QSUB_TEMPLATE_HEADER
            template = QSUB_TEMPLATE

        with open(os.path.join(os.getcwd(), 'qsub_coding_fracs.sh'), 'w') as outfile:
            if opt.user is not None:
                login = opt.user
            else:
                login = os.getlogin()

            num_cores = 2
            max_run = max( [int(i.split('-')[2]) for i in org_out_dirs] )
            outfile.write(header)
            outfile.write(template % (login, 
                                      max_run, ##int(opt.num_runs), 
                                      login,
                                      num_cores, 
                                      #opt.cache_dir, 
                                      opt.features,
                                      opt.organism,
                                      "%s-out-${BATCHNUM}" % (opt.organism),
                                      "%s-out-${BATCHNUM}" % (opt.organism)))
            

if __name__ == '__main__':
    main()


