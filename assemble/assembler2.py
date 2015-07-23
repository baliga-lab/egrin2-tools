#!/usr/bin/env python
import argparse
import logging

from assemble.makeCorems import CoremMaker
import merge
import sqlite3
import pandas as pd


DESCRIPTION = """assemble2.py - prepare cluster runs"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None # "assembler.log"

class SqliteDB:
    """database interface to sqlite"""
    def __init__(self, conn):
        self.conn = conn
        self.conn.execute('create table if not exists row_row (keyrow_id int, subrow_id int, counts int, weight decimal, backbone_pval decimal)')
        self.conn.execute('create table if not exists corems (density decimal, corem_num int, weighted_density decimal)')
        self.conn.execute('create table if not exists corem_rows (corem_id int, row_id int)')
        self.conn.execute('create table if not exists corem_cols (corem_id int, col_id int, pval decimal)')
        self.conn.execute('create table if not exists corem_edges (corem_id int, row1_id int, row2_id int)')

    def get_row_maps(self):
        row2id = {}
        id2row = {}
        cursor = self.conn.cursor()
        try:
            cursor.execute('select rowid,name from rows')
            for rowid, name in cursor.fetchall():
                row2id[name] = rowid
                id2row[rowid] = name
            return row2id, id2row
        finally:
            cursor.close()
    
    def get_column_maps(self):
        col2id = {}
        id2col = {}
        cursor = self.conn.cursor()
        try:
            cursor.execute('select rowid,name from columns')
            for rowid, name in cursor.fetchall():
                col2id[name] = rowid
                id2col[rowid] = name
            return col2id, id2col
        finally:
            cursor.close()

    def num_row_co_occurence(self, rowname, row2id, id2row):
        """for a given row, return the number of co-occurences with genes in all
        biclusters. The result is a Series using the gene names as indexes and
        counts as values"""
        row_pk = row2id[rowname]
        cursor = self.conn.cursor()
        try:
            cursor.execute('select row_id, count(cluster_id) from bicluster_rows where cluster_id in (select cluster_id from bicluster_rows where row_id=?) group by row_id', [row_pk])
            result = {id2row[row_pk]: count for row_pk, count in cursor.fetchall()}
            return pd.Series(result)
        finally:
            cursor.close()

    def drop_row_rows(self):
        self.conn.execute('delete from row_row')

    def update_row_row(self, keyrow_pk, subrow_pk, data_counts_norm, backbone_pval):
        self.conn.execute('update row_row set weight=?, backbone_pval=? where keyrow_id=? and subrow_id=?', [data_counts_norm, backbone_pval, keyrow_pk, subrow_pk])

    def insert_row_row(self, rowrow_docs):
        for doc in rowrow_docs:
            conn.execute('insert into row_row (keyrow_id,subrow_id,counts,weight,backbone_pval) values (?,?,?,?,?)', [doc['row_ids'][0], doc['row_ids'][1], doc['counts'], doc['weight'], doc['backbone_pval']])

    def insert_corem(self, corem_docs):
        cursor = conn.cursor()
        try:
            for doc in corem_docs:
                cursor.execute('insert into corems (density,corem_num,weighted_density) values (?,?,?)', [doc['density'], doc['corem_id'], doc['weighted_density']])
                corem_id = cursor.lastrowid
                for row_id in doc['rows']:
                    cursor.execute('insert into corem_rows (corem_id,row_id) values (?,?)',
                                   [corem_id, row_id])
                for col_id in doc['cols']:
                    cursor.execute('insert into corem_cols (corem_id,col_id) values (?,?)',
                                   [corem_id, col_id])
                for edge_str in doc['edges']:
                    row1, row2 = edge_str.split('-')
                    cursor.execute('insert into corem_edges (corem_id,row1_id,row2_id) values (?,?,?)',
                                   [corem_id, int(row1), int(row2)])
        finally:
            cursor.close()

    def corem_sizes(self):
        cursor = conn.cursor()
        try:
            cursor.execute('select distinct count(row_id) as corem_size from corem_rows group by corem_id order by corem_size')
            return [row[0] for row in cursor.fetchall()]
        finally:
            cursor.close()


if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="3 letter organism code")
    parser.add_argument('--ratios', required=True)
    parser.add_argument('--targetdb', required=True)
    parser.add_argument('--targetdir', required=True, help="Storage path for MongoDB and corem data")
    parser.add_argument('result_dbs', nargs='*')

    # can be overridden
    parser.add_argument('--backbone_pval', default=0.05, type=float, help="Significance pvalue for gene-gene backbone. Default = 0.05.")
    parser.add_argument('--cores', default=3, type=int, help="Number local cores to use for corem C++ scripts")
    parser.add_argument('--link_comm_score', default=0, type=int, help="Scoring metric for link communities" )
    parser.add_argument('--link_comm_increment', default=0.1, type=float, help="Height increment for cutting agglomerative clustering of link communities" )
    parser.add_argument('--link_comm_density_score', default=5, type=int,  help="Density score for evaluating link communities")
    parser.add_argument('--corem_size_threshold', default=3, type=int, help="Defines minimum corem size. Default = 3." )
    parser.add_argument('--n_resamples', default=10000, type=int, help="Number resamples to compute for corem condition assignment. Default = 10,000")

    args = parser.parse_args()

    merge.merge(args)
    conn = sqlite3.connect(args.targetdb, 15, isolation_level='DEFERRED')
    try:
        corems = CoremMaker(args.organism, SqliteDB(conn), args.backbone_pval, args.targetdir,
                            args.cores, args.link_comm_score,
                            args.link_comm_increment,
                            args.link_comm_density_score,
                            args.corem_size_threshold,
                            args.n_resamples)
        corems.make_corems()
        conn.commit()
    finally:
        conn.close()
