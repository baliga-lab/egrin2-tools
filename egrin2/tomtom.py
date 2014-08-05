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
"""
import argparse
import bz2
import os
import subprocess

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="tomtom.py - run tomtom on cmonkey results")
    parser.add_argument('--dir', default='.')
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--targetdir', required=True)

    args = parser.parse_args()
    if not os.path.exists(args.targetdir):
        os.mkdir(args.targetdir)
    resultdirs = [entry for entry in os.listdir(args.dir)
                if entry.startswith(args.prefix) and os.path.isdir(entry)]
    basenames = []
    for resultdir in resultdirs:
        basename = os.path.basename(resultdir)
        print "processing", resultdir, '...'
        export_motifs.export_run_motifs_to_meme(os.path.join(resultdir, 'cmonkey_run.db'),
                                                args.targetdir, basename,
                                                max_residual=None, max_evalue=None)
        basenames.append(basename)

    run_tomtom(args.targetdir,
               os.path.join(args.targetdir, '%s.meme' % basenames[0]),
               os.path.join(args.targetdir, '%s.meme' % basenames[1]))
