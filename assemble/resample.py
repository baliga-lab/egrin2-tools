#!/usr/bin/env python3

"""
Resample ratios. Add them to col_resample collection in an existing egrin2 MongoDB

Example:

python assembler.py --organism mtu --ratios ./20141130.MTB.all.ratios.csv.gz --targetdir ./ --ncbi_code 83332 --ensembledir ./    --n_resamples 0
"""

import random
import argparse
import os
import itertools
import math
import logging

from pymongo import MongoClient
import numpy as np
import pandas as pd
import sqlite3
import json

import assemble.util as util


DESCRIPTION = """resample.py - prepare brute force random resamples"""


def rsd(vals):
    return abs(np.std(vals) / np.mean(vals))

def make_col_id_rsd_groups(df, n_rows, n_resamples):
    # for speed, setup a large data frame to select from
    df_gb = df.dropna().groupby("col_id")
    sel_df = df_gb.aggregate(lambda g: g.tolist())
    df_rsd = pd.concat([sel_df.apply(lambda x: (rsd(random.sample(x[0], n_rows)),
                                                rsd(random.sample(x[1], n_rows))),
                                     axis=1, broadcast=True) for i in range(0, n_resamples)])
    return df_rsd.groupby(df_rsd.index)

class MongoDB:
    def __init__(self, dbclient):
        self.dbclient = dbclient

    def insert_col_resample(self, doc):
        """lowest_raw and lowest_standardized are lists of decimal values
        """
        self.dbclient['col_resample'].insert(doc)

    def update_col_resample(self, n_rows, col_num, n_resamples, lowest_raw, lowest_standardized):
        self.dbclient["col_resample"].update({"n_rows": n_rows, "col_id": col_num},
                                    {"$set": {"resamples": n_resamples,
                                              "lowest_raw": lowest_raw,
                                              "lowest_standardized": lowest_standardized}})

    def find_col_resamples(self, n_rows, cols):
        return self.dbclient["col_resample"].find({ "n_rows": n_rows, "col_id": {"$in": cols}})

    def find_gene_expressions(self, column_nums):
        return pd.DataFrame(list(self.dbclient.gene_expression.find({"col_id": {"$in": column_nums}},
                                                                    {"_id": 0, "col_id": 1,
                                                                     "raw_expression": 1,
                                                                    "standardized_expression": 1})))

    def get_col_nums(self):
        return pd.DataFrame(list(self.dbclient["col_info"].find({}, {"col_id": 1}))).col_id.tolist()


class SqliteDB:
    def __init__(self, conn):
        self.conn = conn
        self.__create_tables()

    def __create_tables(self):
        """we store arrays of decimals into the lowest_raw_exps and
        lowest_std_exps columns by serializing them with the json.dumps/json.loads functions"""
        self.conn.execute('create table if not exists col_resamples (col_id int, nrows int, nresamples int, lowest_raw_exps text, lowest_std_exps text)')

    def insert_col_resample(self, doc):
        self.conn.execute('insert into col_resamples (col_id,nrows,nresamples,lowest_raw_exps,lowest_std_exps) values (?,?,?,?,?)',
                          [doc['col_id'], doc['n_rows'], doc['resamples'],
                           json.dumps(doc['lowest_raw']), json.dumps(doc['lowest_standardized'])])

    def update_col_resample(self, nrows, col_pk, nresamples, lowest_raw, lowest_std):
        self.conn.execute('update col_resamples set nresamples=?,lowest_raw=?,lowest_std=? where col_id=? and nrows=?',
                          [nresamples, json.dumps(lowest_raw), json.dumps(lowest_std), col_pk, nrows])

    def find_col_resamples(self, nrows, col_ids):
        cursor = self.conn.cursor()
        try:
            in_clause = "(%s)" % ",".join(map(str, col_ids))
            query = 'select rowid,col_id,nrows,nresamples,lowest_raw_exps,lowest_std_exps from col_resamples where nrows=? and col_id in ' + in_clause
            cursor.execute(query, [nrows])
            result = []
            for rowid, col_pk, nrows, nresamples, lowest_raw, lowest_std in cursor.fetchall():
                result.append({'_id': rowid, 'col_id': col_pk, 'n_rows': nrows,
                               'resamples': resamples,
                               'lowest_raw': json.loads(lowest_raw),
                               'lowest_standardized': json.loads(lowest_std)})
        finally:
            cursor.close()
        return result

    def find_gene_expressions(self, column_pks):
        cursor = self.conn.cursor()
        try:
            in_list = '(%s)' % ','.join(map(str, column_pks))
            cursor.execute('select col_id,value,std_value from expr_values where col_id in ' + in_list)
            return pd.DataFrame([{'col_id':  col_id, 'raw_expression': value, 'standardized_expression': std_value}
                                 for col_id, value, std_value in cursor.fetchall()])
        finally:
            cursor.close()

    def get_col_nums(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute('select rowid from columns')
            return [row[0] for row in cursor.fetchall()]
        finally:
            cursor.close()


def __choose_n(dbclient, col, vals, n, add, n_rows, n_resamples, old_records, keepP):
    raw = vals.loc[:, "raw_expression"].copy()
    raw.sort_values()
    standardized = vals.loc[:, "standardized_expression"].copy()
    standardized.sort_values()

    if add:
        d = {
            "n_rows": n_rows,
            "col_id": col,
            "resamples": n_resamples,
            "lowest_raw": raw.iloc[0:n].tolist(),
            "lowest_standardized": standardized.iloc[0:n].tolist()
        }
        dbclient.insert_col_resample(d)
    else:
        # update
        resamples = n_resamples + old_records[col]["resamples"]
        n2keep = round(resamples * keepP)
        ran = raw.tolist() + old_records[col]["lowest_raw"]
        ran.sort()
        ran = ran[0: int(n2keep)]
        ras = standardized.tolist() + old_records[col]["lowest_standardized"]
        ras.sort()
        ras = ras[0: int(n2keep)]
        dbclient.update_col_resample(n_rows, col, resamples, ran, ras)


def update_db_col_resample(dbclient, columns, n_rows, n_resamples, keepP,
                           old_records, add=True):
    n2keep = int(round(n_resamples * keepP))
    # do in batches of 100 so memory usage doesn't get too high
    nbins = int(math.ceil(len(columns) / 100.0))
    bins = util.split_list(columns, nbins)

    logging.info("%d bins created", len(bins))
    for index, b in enumerate(bins):
        logging.info("processing bin %d (of %d)", index, nbins)
        df = dbclient.find_gene_expressions(b)
        if df.shape != (0,0):
            logging.info('making rsd on col ids (n_rows = %d)...', n_rows)
            start_time = util.current_millis()
            df_rsd_gb = make_col_id_rsd_groups(df, n_rows, n_resamples)
            elapsed = util.current_millis() - start_time
            logging.info("make_rsd_groups in %d s.", (elapsed / 1000))

            for i in df_rsd_gb.groups.keys():
                __choose_n(dbclient, int(i), df_rsd_gb.get_group(i), n2keep,
                           add=add, n_rows=n_rows,
                           n_resamples=n_resamples, old_records=old_records, keepP=keepP)
        else:
            logging.info("no gene expressions found")


def col_resample_ind(dbclient, n_rows, cols, n_resamples=1000, keepP=0.1):
    """Resample gene expression for a given number of genes in a particular condition using RSD, brute force."""
    logging.info("Adding brute force resample document for gene set size %d", n_rows)

    # see if a record for this gene set size exists
    old_records = {i["col_id"]: i
                   for i in dbclient.find_col_resamples(n_rows, cols)}

    if old_records is not None:
        logging.info("finding records that need to be updated...")
        cols_to_update = [int(i["col_id"]) for i in old_records.values() if int(i["resamples"]) < n_resamples]

    cols_to_add = [i for i in cols if i not in old_records.keys()]

    if len(cols_to_add) == 0 and len(cols_to_pdate) == 0:
        logging.info("Nothing to add")
    else:
        logging.info("%d records to add and %d records to update...", len(cols_to_add), len(cols_to_update))

    if len(cols_to_add) > 0:
        update_db_col_resample(dbclient, cols_to_add, n_rows, n_resamples, keepP,
                               old_records, add=True)

    if len(cols_to_update) > 0:
        n_resamples -= np.min([i["resamples"] for i in old_records.values()])
        if n_resamples > 0:
            update_db_col_resample(dbclient, cols_to_update, n_rows, n_resamples, keepP,
                                   old_records, add=False)


LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None # "resample.log"


if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser( description=DESCRIPTION )
    parser.add_argument('--host', default='localhost', help="Host for MongoDB" )
    parser.add_argument('--db', required=True, help="Database name")
    parser.add_argument('--n_rows', required=True, type=int, help="Gene set size to test")
    parser.add_argument('--n_resamples', default=1000, type=int, help="Number of resamples to compute")
    parser.add_argument('--port', default=27017, help="MongoDB port", type=int )
    parser.add_argument('--keep_p', default=0.1, help="Lowest percent to store", type=int )
    parser.add_argument('--dbengine', default='sqlite', help="Database Engine (sqlite|mongodb)")

    args = parser.parse_args()
    if args.dbengine == 'sqlite':
        conn = sqlite3.connect(args.db)
        dbclient = SqliteDB(conn)
    elif args.dbengine == 'mongodb':
        client = MongoClient(host=args.host, port=args.port)
        dbclient = MongoDB(client[args.db])
    else:
        raise Exception('please specify a database engine')

    cols = dbclient.get_col_nums()
    logging.info("Starting resample for n_rows = %d", args.n_rows)
    col_resample_ind(dbclient, n_rows=args.n_rows, cols=cols, n_resamples=args.n_resamples, keepP=args.keep_p)

    if args.dbengine == 'sqlite':
        conn.commit()
        conn.close()
    elif args.dbengine == 'mongodb':
        client.close()
