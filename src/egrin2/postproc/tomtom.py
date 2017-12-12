#!/usr/bin/env python3

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

./egrin2-tools/postproc/tomtom.py --dir <dir> --prefix <prefix> --targetdir <target directory>
"""
import argparse
import bz2
import os
import subprocess

from egrin2.postproc import export_motifs

# These are default parameters for tomtom that we chose for now
EVALUE_CUTOFF = 100
RESID_CUTOFF  = 0.8
DIST_METHOD   = "ed"
Q_THRESHOLD   = 0.5
MIN_OVERLAP   = 4
Q_PSEUDO      = 0.01
T_PSEUDO      = 0.01


MAX_CLUSTER_RESIDUAL = None
MAX_EVALUE = None

# Template for tomtom 4.9.0
QSUB_TEMPLATE = """#!/bin/bash

export LD_LIBRARY_PATH=/tools/lib:/tools/R-3.0.3/lib64/R/lib
export PATH=/tools/bin:${PATH}

#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=8G

PYTHONPATH=egrin2-tools/src python egrin2-tools/src/postproc/tomtom.py --prefix %s --gene %s

tomtom -verbosity 1 -q-thresh %f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f %s %s | bzip2 -c  > %s
"""

# Template for tomtom 4.9.0
QSUB_TEMPLATE_CSH = """#!/bin/csh -f

export LD_LIBRARY_PATH=/tools/lib:/tools/R-3.0.3/lib64/R/lib
setenv PATH /tools/bin:${PATH}
setenv PYTHONPATH egrin2-tools/src

##$ -m be
#$ -q baliga
#$ -j y
#$ -P Bal_%s
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=8G

python egrin2-tools/src/postproc/tomtom.py --prefix %s --gene %s

tomtom -verbosity 1 -q-thresh %f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f %s %s | bzip2 -c > %s
"""

QSUB_SCRIPT_BASH = """#!/bin/bash
for f in %s/*-tomtom.sh; do
  echo "qsub $f"
  qsub $f
done
"""

QSUB_SCRIPT_CSH = """#!/bin/csh
foreach f (`ls %s/*-tomtom.sh`)
echo "qsub $f"
qsub $f
end
"""

QSUB_CONTROL_SCRIPT = QSUB_SCRIPT_BASH

def emit_tomtom_script(targetdir, filepath, prefix, gene, login, q_thresh=Q_THRESHOLD,
                       dist_method=DIST_METHOD,
                       min_overlap=MIN_OVERLAP, q_pseudo=Q_PSEUDO, t_pseudo=T_PSEUDO):
    num_cores = 1
    with open(os.path.join(targetdir, '%s-tomtom.sh' % gene), 'w') as outfile:
        outfile.write(QSUB_TEMPLATE % (login, login, num_cores, prefix, gene, q_thresh,
                                       dist_method, min_overlap, q_pseudo, t_pseudo,
                                       filepath, filepath, '%s-tomtom.tsv.bz2' % filepath))
    with open('qsub_tomtom.sh', 'w') as outfile:
        outfile.write(QSUB_CONTROL_SCRIPT % targetdir)


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
        print(" ".join(command))
        output = subprocess.check_output(command)
        targetname = os.path.basename(targetfile).rstrip('.meme')
        queryname = os.path.basename(queryfile).rstrip('.meme')
        with open(os.path.join(targetdir, '%s_vs_%s.tsv' % (queryname, targetname)), 'w') as outfile:
            outfile.write(output)
    except:
        raise


def main():
    parser = argparse.ArgumentParser(description="tomtom.py - run tomtom on cmonkey results")
    parser.add_argument('--prefix', required=True, help='a common prefix of the result directories')
    parser.add_argument('--csh', action='store_true')
    parser.add_argument('--user', default=None, help='username for qsub')
    parser.add_argument('--gene', default=None,
                        help='run tomtom only on motifs from clusters containing this gene')

    args = parser.parse_args()
    if args.user is not None:
        login = args.user
    else:
        None

    if args.csh:
        QSUB_TEMPLATE = QSUB_TEMPLATE_CSH
        QSUB_CONTROL_SCRIPT = QSUB_SCRIPT_CSH
    else:
        QSUB_CONTROL_SCRIPT = QSUB_SCRIPT_BASH

    target_dir = 'tomtom_out'
    try:
        os.mkdir(target_dir)
    except:
        None

    if args.gene is not None:
        export_motifs.make_meme_files('.', args.prefix, target_dir, args.gene)
    else:
        genes, dbpaths, max_iteration = export_motifs.get_all_genes('.', args.prefix)
        for gene in genes:
            emit_tomtom_script(target_dir, os.path.join(target_dir, '%s.meme' % gene),
                               args.prefix, gene, login)


if __name__ == '__main__':
    main()
