#!/usr/bin/env python3
import argparse
import sqlite3
import glob
import os

import numpy as np
import pandas

"""
This module reads FIMO information from the postprocessing step and
adds it to the merged results database
"""
def process_fimo_file(conn, dirname, f):
    cursor = conn.cursor()
    try:
        fname = os.path.basename(f)
        cluster = int(fname.replace('fimo-out-', '').replace('.bz2', ''))
        # 1. determine the database id of this cluster
        cursor.execute('select bc.rowid from biclusters bc join ensemble_runs er on er.rowid=bc.run_id where dirname=? and cluster_num=?', [dirname, cluster])
        cluster_id = cursor.fetchone()[0]

        fimo_df = pandas.read_csv(f, sep="\t", compression = "bz2")
        fimo_df.rename(columns={
            'matched sequence': 'matched_sequence',
            'sequence name': 'scaffoldId',
            "#pattern name": "motif_num"
        }, inplace=True)
        trans_d = {}
        for i in np.unique(fimo_df.scaffoldId):
            refseq = "_".join(i.split('.')[-2].split("_")[:: -1][0: 2][:: -1])
            cursor.execute('select mo_scaffold_id from genome where refseq=?', [refseq])
            trans_d[i] = cursor.fetchone()[0]  # scaffold id

        # Real processing happens here
        # translate from reqseqs to MicrobesOnline scaffold ids and store in
        # fimo_df's scaffold id field
        trans_v = [trans_d[i] for i in fimo_df.scaffoldId.values]
        fimo_df.scaffoldId = trans_v
        fimo_df["cluster_id"] = cluster_id

        # only keep specific columns
        fimo_df = fimo_df.loc[:, ['scaffoldId', 'start', 'stop', 'strand', 'score',
                                  'p-value', 'in_coding_rgn', 'cluster_id', 'motif_num']]

        # only store hits with p-value lte 1e-5
        fimo_df = fimo_df.loc[fimo_df["p-value"] <= 1e-5, ]
        data = []
        for index, r in fimo_df.iterrows():
            in_coding_rgn = None
            cursor.execute('select rowid from motif_infos where cluster_id=? and motif_num=?',
                           [cluster_id, r['motif_num']])
            motif_id = cursor.fetchone()[0]
            if not np.isnan(r['in_coding_rgn']):
                in_coding_rgn = r['in_coding_rgn']
            data.append([motif_id, r['strand'], r['start'], r['stop'], r['score'], r['p-value'],
                         in_coding_rgn, r['scaffoldId']])
        conn.executemany('insert into fimo (motif_info_id,strand,start,stop,score,pvalue,in_coding_rgn,mo_scaffold_id) values (?,?,?,?,?,?,?,?)', data)

        # The original inserted fimo_small, which is a subset of fimo that only contains the
        # entries whose motif_info entry also contains a GRE. But this is a redundancy
        # and SQL can find these entries easily using a merge
        conn.commit()
    except:
        print("no data, skip")
    finally:
        cursor.close()


def store_fimo(conn, args):
    prefix = '%s-out-' % args.organism
    cmout_dirs = sorted(glob.glob(os.path.join(args.ensembledir, "%s???" % prefix)))
    for cmout in cmout_dirs:
        print("processing ", cmout)
        dirname = cmout.split('/')[-1]
        fimo_outs = sorted(os.listdir(os.path.join(cmout, 'fimo-outs')))
        for f in fimo_outs:
            process_fimo_file(conn, dirname, os.path.join(cmout, 'fimo-outs', f))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="fimo.py - add fimo information")
    parser.add_argument('organism')
    parser.add_argument('targetdb')
    parser.add_argument('--ensembledir', default='.', help="Path to ensemble runs. Default: cwd")
    args = parser.parse_args()
    conn = sqlite3.connect(args.targetdb)
    store_fimo(conn, args)
