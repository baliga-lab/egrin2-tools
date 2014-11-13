#!/usr/bin/env python

"""genebygene.py - calculate counts for gene-gene relationships

This script scans a directory of ensemble runs and extracts cluster-gene relationships.
It then generates a structure that is stored as a file in JSON format
"""
import argparse
import sqlite3
import itertools
import os
from collections import defaultdict
import json
import pickle
import numpy as np

def process_dir(resultdir, ggcounts):
    print resultdir
    cluster_genes = defaultdict(set)
    dbfile = os.path.join(resultdir, 'cmonkey_run.db')
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from row_members')
    iteration = cursor.fetchone()[0]
    cursor.execute('select cluster, name from row_members rm join row_names rn on rm.order_num = rn.order_num where iteration=?', [iteration])
    for cluster, name in cursor.fetchall():
        cluster_genes[cluster].add(name)
    for cluster, genes in cluster_genes.items():
        genes = sorted(genes)
        for pair in itertools.product(genes, genes):
            if pair[0] != pair[1]:
                ggcounts[pair] += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="genebygene.py - generate gene-by-gene counts")
    parser.add_argument('--dir', default='.', help="directory holding the ensemble run results")
    parser.add_argument('--prefix', required=True, help='a common prefix of the result directories')
    parser.add_argument('--outfile', required=True, help='output matrix file')
    args = parser.parse_args()

    inpath = args.dir
    prefix = args.prefix
    resultdirs = map(lambda s: os.path.join(inpath, s),
                     sorted([entry for entry in os.listdir(inpath)
                             if entry.startswith(prefix) and
                             os.path.isdir(os.path.join(inpath, entry))]))

    ggcounts = defaultdict(int)
    conn = sqlite3.connect(os.path.join(resultdirs[0], 'cmonkey_run.db'))
    cursor = conn.cursor()
    cursor.execute('select name from row_names')
    genenames = [row[0] for row in cursor.fetchall()]
    cursor.close()
    conn.close()

    for index, resultdir in enumerate(resultdirs):
        process_dir(resultdir, ggcounts)

    #with open('gbg.pkl', 'w') as outfile:
    #    pickle.dump(ggcounts, outfile)

    genenames = sorted(genenames)
    num_genes = len(genenames)
    mat = np.zeros((num_genes, num_genes))
    for i, rgene in enumerate(genenames):
        for j, cgene in enumerate(genenames):
            if (rgene, cgene) in ggcounts:
                mat[i, j] = ggcounts[(rgene, cgene)]
            elif (cgene, rgene) in ggcounts:
                mat[i, j] = float(ggcounts[(cgene, rgene)])

    # normalize the matrix so the values sum to 1
    row_sums = mat.sum(axis=1)
    mat = (mat.T / row_sums).T
    
    with open(args.outfile, 'w') as outfile:
        titlerow = '\t'.join(genenames)
        outfile.write(titlerow + '\n')
        for i, rgene in enumerate(genenames):
            row = [rgene]
            row.extend([('%.5f' % value) for value in mat[i]])
            outfile.write('\t'.join(row) + '\n')
