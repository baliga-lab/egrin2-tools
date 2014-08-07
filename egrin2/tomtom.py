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
import sqlite3


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


# Template for MEME file header

MEME_FILE_HEADER = """MEME version 3.0

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
A %.3f C %.3f G %.3f T %.3f
"""


def write_pssm(outfile, cursor, run_name, cluster, motif_info_id,
               motif_num, evalue, num_sites):
    """writes a single PSSM to the given file"""
    motif_name = '%s_%03d_%02d' % (run_name, cluster, motif_num)
    outfile.write('\nMOTIF %s\n' % motif_name)
    outfile.write('BL   MOTIF %s width=0 seqs=0\n' % motif_name)

    cursor.execute('select a,c,g,t from motif_pssm_rows where motif_info_id=? order by row',
                   [motif_info_id])
    pssm_rows = [(a, c, g, t) for a, c, g, t in cursor.fetchall()]
    outfile.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
    for a, c, g, t in pssm_rows:
        outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (a, c, g, t))


def make_meme_file(dbpaths, maxiter, targetdir, gene,
                   max_residual=None, max_evalue=None,
                   a_perc=0.284, c_perc=0.216, g_perc=0.216, t_perc=0.284):
    """Creates a meme file for a specific gene. Returns the number of
    PSSMs that were written"""
    num_written = 0
    targetpath = os.path.join(targetdir, '%s.meme' % gene)
    with open(targetpath, 'w') as outfile:
        outfile.write(MEME_FILE_HEADER % (a_perc, c_perc, g_perc, t_perc))
        # gene -> all runs
        for dbpath in dbpaths:
            conn = sqlite3.connect(dbpath)
            cursor = conn.cursor()
            cursor2 = conn.cursor()
            cursor.execute('select order_num from row_names where name=?', [gene])
            order_num = cursor.fetchone()[0]
            query = """select mi.rowid, mi.cluster, motif_num, evalue,
                       count(mms.pvalue) as num_sites
                       from (select rowid,cluster,motif_num,evalue from motif_infos
                             where iteration=?) mi
                       join (select cluster from row_members where iteration=?
                              and order_num=?) as rm on mi.cluster=rm.cluster
                       join (select cluster, residual from cluster_residuals
                             where iteration=?) cr on rm.cluster=cr.cluster
                       join meme_motif_sites mms on mi.rowid=mms.motif_info_id"""
            params = [maxiter, maxiter, order_num, maxiter]

            # optional residual and evalue filters
            if max_residual is not None or max_evalue is not None:
                query += " where"
            if max_residual is not None:
                query += " residual <= ?"
                params.append(max_residual)
            if max_evalue is not None:
                if max_residual is not None:
                    query += " and "
                query += " evalue <= ?"
                params.append(max_evalue)

            query += " group by mi.rowid"
            cursor.execute(query, params)
            for rowid, cluster, motif_num, evalue, num_sites in cursor.fetchall():
                write_pssm(outfile, cursor2,
                           os.path.basename(os.path.dirname(dbpath)),
                           cluster, rowid, motif_num, evalue, num_sites)
                num_written += 1
            cursor2.close()
            cursor.close()
            conn.close()
    return num_written


def make_meme_files(inpath, prefix, targetdir):
    """create MEME files based on each gene in the ensemble run and writing
    all clusters in all runs that contain the gene"""
    resultdirs = map(lambda s: os.path.join(inpath, s),
                     sorted([entry for entry in os.listdir(inpath)
                             if entry.startswith(prefix) and os.path.isdir(entry)]))
    dbpaths = [os.path.join(resultdir, 'cmonkey_run.db') for resultdir in resultdirs]
    # extract max iteration
    conn = sqlite3.connect(dbpaths[0])
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from row_members')
    max_iteration = cursor.fetchone()[0]
    cursor.execute('select name from row_names order by name')
    genes = [row[0] for row in cursor.fetchall()]
    conn.close()
    for gene in genes:
        print "processing gene '%s'..." % gene,
        num_written = make_meme_file(dbpaths, max_iteration, targetdir, gene)
        print "%d motifs written." % num_written


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="tomtom.py - run tomtom on cmonkey results")
    parser.add_argument('--dir', default='.')
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--targetdir', required=True)

    args = parser.parse_args()
    if not os.path.exists(args.targetdir):
        os.mkdir(args.targetdir)
    make_meme_files(args.dir, args.prefix, args.targetdir)

    """
    run_tomtom(args.targetdir,
               os.path.join(args.targetdir, '%s.meme' % basenames[0]),
               os.path.join(args.targetdir, '%s.meme' % basenames[1]))"""
