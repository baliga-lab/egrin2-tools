#!/usr/bin/env python

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

from query.egrin2_query import *

DESCRIPTION = """resample.py - prepare brute force random resamples"""


def rsd(vals):
    return abs(np.std(vals) / np.mean(vals))


def resample(row_vals, n_rows):
    """filter out nan values!!!!"""
    return rsd(random.sample(row_vals.dropna(), n_rows))

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


def __choose_n(dbclient, col, vals, n, add, n_rows, n_resamples, old_records, keepP):
    raw = vals.loc[:, "raw_expression"].copy()
    raw.sort()
    standardized = vals.loc[:, "standardized_expression"].copy()
    standardized.sort()

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


def __split_list(alist, wanted_parts=1):
    length = len(alist)
    return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
            for i in range(wanted_parts)]

def col_resample_ind(dbclient, n_rows, cols, n_resamples=1000, keepP=0.1):
    """Resample gene expression for a given number of genes in a particular condition using RSD, brute force."""
    logging.info("Adding brute force resample document for gene set size %d", n_rows)

    # see if a record for this gene set size exists
    old_records = {i["col_id"]: i
                   for i in dbclient.find_col_resamples(n_rows, cols)}

    if old_records is not None:
        logging.info("finding records that need to be updated...")
        toUpdate = [int(i["col_id"])
                    for i in old_records.values() if int(i["resamples"]) < n_resamples]

    toAdd = [i for i in cols if i not in old_records.keys()]

    if len(toAdd) == 0 and len(toUpdate) == 0:
        logging.info("Nothing to add")
        return None
    else:
        logging.info("%d records to add and %d records to update...", len(toAdd), len(toUpdate))


    n2keep = int(round(n_resamples * keepP))

    # toAdd
    if len(toAdd) > 0:
        logging.info("Computing resamples for new entries")

        # do in batches of 100 so memory usage doesn't get too high
        nbins = int(math.ceil(len(toAdd) / 100.0))
        bins = __split_list(toAdd, nbins)
        logging.info("%d bins created", len(bins))
        for index, b in enumerate(bins):
            logging.info("processing bin %d", index)
            df = dbclient.find_gene_expressions(b)
            if df.shape != (0,0):
                df_gb = df.groupby("col_id")
                logging.info('making rsd on col ids (n_rows = %d)...', n_rows)
                df_rsd = pd.concat([df_gb.aggregate(resample, n_rows) for i in xrange(0, n_resamples)])
                df_rsd_gb = df_rsd.groupby(df_rsd.index)

                logging.info("Adding new entries...")
                for i in df_rsd_gb.groups.keys():
                    __choose_n(dbclient, int(i), df_rsd_gb.get_group(i), n2keep, True, n_rows,
                               n_resamples, old_records, keepP)

    # toUpdate
    if len(toUpdate) > 0:
        logging.info("Computing resamples for updated entries")

        # do in batches of 500 so memory usage doesn't get too high
        nbins = int(math.ceil(len(toUpdate) / 100.0))
        bins = split_list(toUpdate, nbins)
        for b in bins:
            df = dbclient.find_gene_expressions(self, column_nums)

            if df.shape != (0, 0):
                df_gb = df.groupby("col_id")
                resamples = n_resamples - np.min([i["resamples"] for i in old_records.values()])

                if resamples > 0:
                    df_rsd = pd.concat([df_gb.aggregate(resample, n_rows) for i in xrange(0, resamples)])
                    df_rsd_gb = df_rsd.groupby(df_rsd.index)

                    logging.info("Updating entries")
                    for i in df_rsd_gb.groups.keys():
                        __choose_n(dbclient, int(i), df_rsd_gb.get_group(i), n2keep, False, n_rows,
                                   resamples, old_records, keepP)
    return None


LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None # "resample.log"


if __name__ == '__main__':
    import argparse
    import os
    import itertools

    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser( description=DESCRIPTION )
    parser.add_argument('--host', required=True, type=str, help="Host for MongoDB" )
    parser.add_argument('--db', required=True, type=str, help="Database name")
    parser.add_argument('--n_rows', required=True, type=int, help="Gene set size to test")
    parser.add_argument('--n_resamples', default=1000, type=int, help="Number of resamples to compute")
    parser.add_argument('--port', default=27017, help="MongoDB port", type=int )
    parser.add_argument('--keep_p', default=0.1, help="Lowest percent to store", type=int )
    parser.add_argument('--cols', default=None, help="Columns (experiments) for resampling. Should be path to tab-delimited file containing column names that map to egrin2_col_names in MongoDB database, eg experiment names.", type=str )

    args = parser.parse_args()
    client = MongoClient(host=args.host, port=args.port)
    dbclient = MongoDB(client[args.db])

    if args.cols is None:
        cols = dbclient.get_col_nums()
    else:
        # not supported yet
        logging.error("Not supported yet")
        cols = None

    logging.info("Starting resample for n_rows = %d", args.n_rows)
    col_resample_ind(dbclient, n_rows=args.n_rows, cols=cols, n_resamples=args.n_resamples, keepP=args.keep_p)
    client.close()
