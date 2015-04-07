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


def choose_n(col, vals, n, add, client, db, n_rows, n_resamples, old_records, keepP):
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
        client[db]["col_resample"].insert(d)
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
        client[db]["col_resample"].update({"n_rows": n_rows, "col_id": col},
                                          {"$set": {"resamples": resamples, "lowest_raw": ran, "lowest_standardized": ras}})

def colResampleInd(host, db, n_rows, cols, n_resamples=1000, keepP=0.1, port=27017):
    """Resample gene expression for a given number of genes in a particular condition using RSD, brute force."""

    logging.info("Adding brute force resample document for gene set size %d", n_rows)

    # make connection
    client = MongoClient(host=host, port=port)

    # see if a record for this gene set size exists
    old_records = {i["col_id"]: i
                   for i in client[db]["col_resample"].find({ "n_rows": n_rows, "col_id": {"$in": cols}})}

    if old_records is not None:
        toUpdate = [int(i["col_id"])
                    for i in old_records.values() if int(i["resamples"]) < n_resamples]

    toAdd = [i for i in cols if i not in old_records.keys()]

    if len(toAdd) == 0 and len(toUpdate) == 0:
        logging.info("Nothing to add")
        client.close()
        return None

    n2keep = int(round(n_resamples * keepP))

    def split_list(alist, wanted_parts=1):
        length = len(alist)
        return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
                for i in range(wanted_parts)]

    # toAdd
    if len(toAdd) > 0:
        logging.info("Computing resamples for new MongoDB documents")

        # do in batches of 100 so memory usage doesn't get too high
        nbins = int(math.ceil(len(toAdd) / 100.0))
        bins = split_list(toAdd, nbins)
        for b in bins:
            df = pd.DataFrame(list(client[db].gene_expression.find({"col_id": {"$in": b}},
                                                                   { "col_id": 1, "raw_expression": 1, "standardized_expression": 1})))
            if df.shape != (0,0):
                df = df.groupby("col_id")
                df_rsd = pd.concat([df.aggregate(resample, n_rows) for i in range(0, n_resamples)])
                df_rsd = df_rsd.groupby(df_rsd.index)

                logging.info("Adding new documents to MongoDB")
                tmp = [choose_n(int(i), df_rsd.get_group(i), n2keep, True, client, db, n_rows,
                                n_resamples, old_records, keepP) for i in df_rsd.groups.keys()]

    # toUpdate
    if len(toUpdate) > 0:
        logging.info("Computing resamples for updated MongoDB documents")

        # do in batches of 500 so memory usage doesn't get too high
        nbins = int(math.ceil(len(toUpdate) / 100.0))
        bins = split_list(toUpdate, nbins)
        for b in bins:
            df = pd.DataFrame(list( client[db].gene_expression.find({ "col_id": {"$in": b}},
                                                                    {"col_id":1, "raw_expression":1, "standardized_expression": 1})))

            if df.shape != (0, 0):
                df = df.groupby("col_id")
                resamples = n_resamples - np.min([i["resamples"] for i in old_records.values()])

                if resamples > 0:
                    df_rsd = pd.concat([df.aggregate(resample, n_rows) for i in range(0, resamples)])
                    df_rsd = df_rsd.groupby(df_rsd.index)

                    logging.info("Updating MongoDB documents")
                    tmp = [choose_n(int(i), df_rsd.get_group(i), n2keep, False, client, db, n_rows,
                                    resamples, old_records, keepP) for i in df_rsd.groups.keys()]
    client.close()

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
    parser.add_argument( '--host', required=True, type=str, help="Host for MongoDB" )
    parser.add_argument('--db', required=True, type=str, help="Database name")
    parser.add_argument('--n_rows', required=True, type=int, help="Gene set size to test")
    parser.add_argument('--n_resamples', default=1000, type=int, help="Number of resamples to compute")
    parser.add_argument('--port', default=27017, help="MongoDB port", type=int )
    parser.add_argument('--keep_p', default=0.1, help="Lowest percent to store", type=int )
    parser.add_argument('--cols', default=None, help="Columns (experiments) for resampling. Should be path to tab-delimited file containing column names that map to egrin2_col_names in MongoDB database, eg experiment names.", type=str )

    args = parser.parse_args()
    client = MongoClient(host=args.host, port=args.port)

    if args.cols is None:
        cols = pd.DataFrame(list(client[args.db]["col_info"].find({}, {"col_id": 1}))).col_id.tolist()
    else:
        # not supported yet
        logging.error("Not supported yet")
        cols = None

    logging.info("Starting resample for n_rows = %d", args.n_rows)

    colResampleInd(host=args.host, db=args.db, n_rows=args.n_rows, cols=cols, n_resamples=args.n_resamples, keepP=args.keep_p, port=args.port)
    client.close()
