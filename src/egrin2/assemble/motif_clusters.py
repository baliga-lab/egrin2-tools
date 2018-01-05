#!/usr/bin/env python3
"""motif_clusters.py - add motif clustering information to the result database

This takes the datafile from the motif clustering tool that was run in the
post processing step and determines the GREs which get stored into the
GRE field of the motif_infos table
"""
import os
import argparse
import sqlite3

DEFAULT_MOTIF_OUT_FILENAME = "out.mot_metaclustering.txt.I24"

def load_gre_map(ensembledir):
    """read mapping:
    {run name -> {cluster -> { motif: gre_id}}}
    """
    filepath = os.path.join(ensembledir, DEFAULT_MOTIF_OUT_FILENAME)
    gre_id = 1  # this is simply a global unique integer
    mots = {}
    with open(filepath, 'r') as f:
        for line in f:
            # only consider motif clusters with > 3 motifs
            if len(line.strip("\n").split("\t")) > 3:
                for motif in line.strip("\n").split("\t"):
                    elements = motif.split("_")
                    # cluster
                    elements[1] = int(elements[1])
                    # motif_num
                    elements[2] = int(elements[2])
                    # elements[0] = run_name
                    if elements[0] in mots.keys():
                        if elements[1] in mots[elements[0]].keys():
                            mots[elements[0]][elements[1]][elements[2]] = gre_id
                        else:
                            mots[elements[0]][elements[1]] = {}
                            mots[elements[0]][elements[1]][elements[2]] = gre_id
                    else:
                        mots[elements[0]] = {}
                        mots[elements[0]][elements[1]] = {}
                        mots[elements[0]][elements[1]][elements[2]] = gre_id
                gre_id += 1
    return mots


def store_gres(conn, motif2gre):
    cursor = conn.cursor()
    update_query = "update motif_infos set gre_id=? where cluster_id=? and motif_num=?"
    try:
        rundirs = sorted(motif2gre.keys())
        for rundir in rundirs:
            clusters = motif2gre[rundir]
            cursor.execute('select rowid from ensemble_runs where dirname=?', [rundir])
            run_id = cursor.fetchone()[0]
            print("run: %s, id: %d, %d clusters..." % (rundir, run_id, len(clusters)), end="")
            cursor.execute('select rowid,cluster_num from biclusters where run_id=?',
                           [run_id])
            cluster_ids = {cluster_num: cluster_id
                           for cluster_id, cluster_num in cursor.fetchall()}
            for cluster_num, motifs in clusters.items():
                for motif_num, gre_id in motifs.items():
                    cursor.execute(update_query, [gre_id, cluster_ids[cluster_num], motif_num])
            conn.commit()
            print("done.")
    finally:
        cursor.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="motif_clusters.py - add mcl information")
    parser.add_argument('ensembledir', help="Path to ensemble runs.")
    parser.add_argument('targetdb', help="path to SQLite database of the assembled results")
    args = parser.parse_args()
    conn = sqlite3.connect(args.targetdb)
    motif2gre = load_gre_map(args.ensembledir)
    #print(motif2gre)
    store_gres(conn, motif2gre)
