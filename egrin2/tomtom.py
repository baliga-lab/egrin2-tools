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


def run_tomtom(targetdir, targetfile, queryfile, q_thresh=Q_THRESHOLD, dist_method=DIST_METHOD,
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
                                                args.targetdir, basename)
    
        basenames.append(basename)

    run_tomtom(args.targetdir,
               os.path.join(args.targetdir, '%s.meme' % basenames[0]),
               os.path.join(args.targetdir, '%s.meme' % basenames[1]))
