#!/usr/bin/env python

import argparse
import bz2
import os
import subprocess

import export_motifs

EVALUE_CUTOFF = 100
RESID_CUTOFF  = 0.8
DIST_METHOD   = "ed"
Q_THRESHOLD   = 0.5
MIN_OVERLAP   = 4
Q_PSEUDO      = 0
T_PSEUDO      = 0


def run_tomtom(targetfile, queryfile, q_thresh=Q_THRESHOLD, dist_method=DIST_METHOD,
               min_overlap=MIN_OVERLAP, q_pseudo=Q_PSEUDO, t_pseudo=T_PSEUDO):
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
        output = subprocess.check_output(command)
        print output
    except:
        print " ".join(command)
        raise

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="tomtom.py - run tomtom on cmonkey results")
    parser.add_argument('--resultdb', required=True)
    parser.add_argument('--targetdir', required=True)

    args = parser.parse_args()
    if not os.path.exists(args.targetdir):
        os.mkdir(args.targetdir)

    export_motifs.export_run_motifs_to_meme(args.resultdb, args.targetdir, 'myrun')

    #run_tomtom(os.path.join(args.targetdir, 'cmresults-postproc-15.meme'),
    #           os.path.join(args.targetdir, 'cmresults-postproc-142.meme'))
