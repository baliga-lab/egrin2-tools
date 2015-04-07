#!/usr/bin/env

import argparse
import logging
import pymongo
from collections import namedtuple

import query.egrin2_query as e2q
import pandas as pd


DESCRIPTION = """assemble_finish.py - finish and dump mongodb"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None


def compute_and_write_col(db, corem, cond_ids, dbname, host, port, n_resamples=1000):
    logging.info("Adding conditions for corem %d", corem['corem_id'])
    pvals = e2q.col_resample_pval(corem['rows'], "row_id", cond_ids, "col_id",
                                  n_resamples, host, port, dbname, keepP=0.05)

    pvals["col_id"] = pvals.index
    d = pvals.to_dict('records')
    print "corem %s -> %s" % (str(corem['_id']), str(d))
    # TODO: enable this line when everything works
    #self.db.corem.update({"_id": corem._id}, {"$set": {"cols": d}})


def finish_corems(db, dbname, host, port):
    """Finish adding corem info (cols) after resampling. Assumes corem docs already exist"""
    cond_ids = [entry['col_id'] for entry in db["col_info"].find({}, {"_id": 0, "col_id": 1})]

    for corem in db["corem"].find({}, {'_id': 1, "rows": 1, "corem_id": 1}):
        compute_and_write_col(db, corem, cond_ids, dbname, host, port)


if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('db', help="MongoDB database name")
    parser.add_argument('--host', default='localhost', required=False, help="MongoDB database host")
    parser.add_argument('--port', default=27017, type=int, required=False, help="MongoDB database host")
    args = parser.parse_args()
    print "Connecting to database: ", args.db
    client = pymongo.MongoClient(host=args.host, port=args.port)
    try:
        finish_corems(client[args.db], args.db, args.host, args.port)
    finally:
        client.close()
