#!/usr/bin/env python

"""tomtom.py - running MEME tomtom on cmonkey-python ensemble runs

This script assumes an ensemble run under a directory structure that 
is as follows:

<dir>
  |-- <prefix>001
             |-- cmonkey_run.db
             |-- ...
  |-- <prefix>002
  |-- ...

It extracts the PSSMs from the clusters in each individual ensemble run
and then runs them against the motifs of enother ensemble run.

The results of tomtom are stored in an appropriately named tab-separated
file in the target directory

Usage:

./tomtom.py --dir <dir> --prefix <prefix> --targetdir <target directory>
"""
import argparse
import bz2
import os
import subprocess
import sqlite3

import export_motifs

# These are default parameters for tomtom that we chose for now
EVALUE_CUTOFF = 100
RESID_CUTOFF  = 0.8
DIST_METHOD   = "ed"
Q_THRESHOLD   = 0.5
MIN_OVERLAP   = 4
Q_PSEUDO      = 0
T_PSEUDO      = 0


MAX_CLUSTER_RESIDUAL = None
MAX_EVALUE = None

# Template for tomtom 4.9.0
QSUB_TEMPLATE = """#!/bin/bash

export PATH=/tools/bin:${PATH}

#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

tomtom -verbosity 1 -q-thresh %f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f %s %s > %s"""

# Template for tomtom 4.9.0
QSUB_TEMPLATE_CSH = """#!/bin/csh

setenv PATH /tools/bin:${PATH}

#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

tomtom -verbosity 1 -q-thresh %f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f %s %s > %s"""

def run_tomtom(targetdir, targetfile, queryfile, q_thresh=Q_THRESHOLD, dist_method=DIST_METHOD,
               min_overlap=MIN_OVERLAP, q_pseudo=Q_PSEUDO, t_pseudo=T_PSEUDO):
    """a wrapper around the tomtom script"""
    command = ['tomtom',
               '-verbosity', '1',
               '-q-thresh', '%.3f' % q_thresh,
               '-dist', dist_method,
               '-min-overlap', '%d' % min_overlap,
               '-text',
               '-query-pseudo', '%.3f' % q_pseudo,
               '-target-pseudo', '%.3f' % t_pseudo,
               '-target', targetfile, '-query', queryfile]
    try:
        print " ".join(command)
        output = subprocess.check_output(command)
        targetname = os.path.basename(targetfile).rstrip('.meme')
        queryname = os.path.basename(queryfile).rstrip('.meme')
        with open(os.path.join(targetdir, '%s_vs_%s.tsv' % (queryname, targetname)), 'w') as outfile:
            outfile.write(output)
    except:
        raise

def emit_tomtom_script(targetdir, filepath, gene, q_thresh=Q_THRESHOLD, dist_method=DIST_METHOD,
               min_overlap=MIN_OVERLAP, q_pseudo=Q_PSEUDO, t_pseudo=T_PSEUDO):
    login = 'mharris'
    num_cores = 1

    with open(os.path.join(targetdir, 'cluster_tomtom-%s.sh' % gene), 'w') as outfile:
        outfile.write(QSUB_TEMPLATE % (login, login, num_cores, q_thresh,
                                       dist_method, min_overlap, q_pseudo, t_pseudo,
                                       filepath, filepath, '%s-tomtom.tsv' % filepath))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="tomtom.py - run tomtom on cmonkey results")
    parser.add_argument('--dir', default='.', help="directory holding the ensemble run results")
    parser.add_argument('--prefix', required=True, help='a common prefix of the result directories')
    parser.add_argument('--targetdir', required=True, help='the directory to store the results')
    parser.add_argument('--csh', action='store_true')

    args = parser.parse_args()

    if args.csh:
      QSUB_TEMPLATE = QSUB_TEMPLATE_CSH

    if not os.path.exists(args.targetdir):
        os.mkdir(args.targetdir)
    genes = export_motifs.make_meme_files(args.dir, args.prefix, args.targetdir)    
    for gene in genes:
        emit_tomtom_script(args.targetdir, os.path.join(args.targetdir, '%s.meme' % gene), gene)
