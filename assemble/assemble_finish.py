#!/usr/bin/env python

import argparse
import logging
import pymongo
import sqlite3
from collections import namedtuple
import itertools

import query.egrin2_query as e2q
import assemble.resample as resample
import pandas as pd


DESCRIPTION = """assemble_finish.py - finish and dump mongodb"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None


class SqliteDB:
    def __init__(self, conn):
        self.conn = conn

    def close(self):
        self.conn.close()

    def get_cond_ids(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute('select rowid from columns')
            return [row[0] for row in cursor.fetchall()]
        finally:
            cursor.close()

    def get_corems(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute('select corem_id, row_id from corem_rows')
            corem_rows = [(corem_id, row_id) for corem_id, row_id in cursor.fetchall()]
            result = [{'_id': corem_id, 'corem_id': corem_id, 'rows': map(lambda x: x[1], row_ids)}
                      for corem_id, row_ids in itertools.groupby(corem_rows, lambda x: x[0])]
            return result
        finally:
            cursor.close()

class MongoDB:
    def __init__(self, dbclient):
        self.dbclient = dbclient

    def close(self):
        self.dbclient.close()

    def get_cond_ids(self):
        return [entry['col_id'] for entry in self.dbclient["col_info"].find({}, {"_id": 0, "col_id": 1})]

    def get_corems(self):
        return [corem for corem in self.dbclient["corem"].find({}, {'_id': 1, "rows": 1, "corem_id": 1})]

    def check_col_resamples(self, col, n_rows, n_resamples):
        return self.dbclient.col_resample.find_one({"n_rows": n_rows, "col_id": col, "resamples": {"$gte": n_resamples}}) is None

    def find_gene_expressions(self, rows, cols):
        return list(self.dbclient.gene_expression.find({"col_id": {"$in": cols}, "row_id": {"$in": rows}},
                                                    { "_id": 0, "col_id": 1, "normalized_expression": 1,
                                                        "standardized_expression": 1}))

    def find_col_resamples(self, rows, cols):
        return list(self.dbclient.col_resample.find({"n_rows": len(rows), "col_id": {"$in": cols}}, {"_id": 0}))

    def get_col2ids_colout(self, cols, col_outtype):
        return e2q.col2id_batch(self.dbclient, cols, input_type="col_id", return_field=col_outtype)

    def update_corem(self, corem, new_cols):
        #print "update_corem()", corem, " conds: ", new_cols
        self.dbclient['corem'].update({"_id": corem['_id']}, {"$set": {"cols": new_cols}})


def __col_resample_pval(dbclient, rows, cols, n_resamples,
                        standardized=True, sig_cutoff=0.05,
                        sort=True, add_override=False, n_jobs=4, keepP=0.1, verbose=False,
                        col_outtype="col_id"):

    def empirical_pval(i, random_rsd, resamples):
        for x in range(0, len(i)):
            val = (float(sum([i.values[0] >= y for y in random_rsd[i.index[0]]])) /
                   len(random_rsd[i.index[0]])) * (float(len(random_rsd[i.index[0]])) /
                                                   resamples[i.index[0]])
            if val >= float(len(random_rsd[i.index[0]])) / resamples[i.index[0]]:
                return round(float(len(random_rsd[i.index[0]])) / resamples[i.index[0]], 2)
            elif val == 0:
                return 1.0 / resamples[i.index[0]]
            else:
                return val

    if len(rows) == 0:
        logging.info("Please provide an appropriately named array of rows")
        return None

    if len(cols) == 0:
        logging.info("Please provide an appropriately named array of cols")
        return None

    # Determine what/how many resamples need to be added to db
    to_add = filter(lambda col: dbclient.check_col_resamples(col, len(rows), n_resamples), cols)

    count = 1
    if len(to_add) > 0:
        if add_override:
            logging.info("I need to perform %d random resample(s) of size %d to compute pvals. Please be patient. This may take a while...", len(to_add), n_resamples)
            tmp = resample.col_resample_ind(dbclient, n_rows=len(rows), cols=to_add, n_resamples=n_resamples, keepP=keepP)
            logging.info("Done adding random resamples.")
        else:
            logging.info("""I would need to perform %d random resample(s) of size %d
to compute pvals. Since this would require significant computational power (and time),
I have only returned results where resample data has been pre-calculated.
Consult resample.py to run these jobs on multiple cores (much faster)
or change 'add_override' flag of this function to 'True' to build the resample now.""",
                         len(to_add), n_resamples)
            cols = [i for i in cols if i not in to_add]

    if verbose:
        logging.info("Calculating pvals...")

    exp_df = pd.DataFrame(dbclient.find_gene_expressions(rows, cols))
    random_rsd = pd.DataFrame(dbclient.find_col_resamples(rows, cols))

    if random_rsd.shape[0] == 0:
        logging.info("Could not find resample DB entry for %d rows in cols %s",
                     len(rows_o), cols)
        return None
    else:
        random_rsd.index = random_rsd["col_id"]

    exp_df_rsd = exp_df.groupby("col_id").aggregate(e2q.rsd)

    if standardized:
        exp_df_rsd = exp_df_rsd.loc[ :, "standardized_expression"]
        resamples = random_rsd.loc[ :, "resamples"].to_dict()
        random_rsd = random_rsd.loc[ :,"lowest_standardized" ].to_dict()
    else:
        exp_df_rsd = exp_df_rsd.loc[ :, "raw_expression"]
        resamples = random_rsd.loc[ :, "resamples"].to_dict()
        random_rsd = random_rsd.loc[ :,"lowest_raw"].to_dict()

    pvals = exp_df_rsd.groupby(level=0).aggregate(empirical_pval, random_rsd, resamples)
    pvals.columns = ["pval"]
    pvals.index = dbclient.get_col2ids_colout(pvals.index.values, col_outtype)

    if sig_cutoff is not None:
        pvals = pvals[pvals <= sig_cutoff]

    if sort:
        pvals.sort()

    pvals = pvals.to_frame()
    pvals.columns = ["pval"]

    if pvals.shape[0] == 0:
        logging.info("No cols pass the significance cutoff of %f", sig_cutoff)

    return pvals


def __compute_and_write_col(dbclient, corem, cond_ids, n_resamples=1000):
    logging.info("Adding conditions for corem %d", corem['corem_id'])
    pvals = __col_resample_pval(dbclient, corem['rows'], cond_ids, n_resamples, keepP=0.05)
    if pvals is not None:
        pvals["col_id"] = pvals.index
        d = pvals.to_dict('records')
        print "corem %s\n" % (str(corem['_id']))
        dbclient.update_corem(corem, d)


def finish_corems(dbclient):
    """Finish adding corem info (cols) after resampling. Assumes corem docs already exist"""
    cond_ids = dbclient.get_cond_ids()
    # read the corems in advance, so we do not hold cursors for too long
    corems = dbclient.get_corems()

    for corem in corems:
        try:
            __compute_and_write_col(dbclient, corem, cond_ids)
        except:            
            logging.exception('ERROR on corem %d %s' % (corem['corem_id'], str(corem['_id'])))
            raise
            

if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('db', help="database name")
    parser.add_argument('--host', default='localhost', help="MongoDB database host")
    parser.add_argument('--port', default=27017, type=int, help="MongoDB database port")
    parser.add_argument('--dbengine', default='sqlite', help="Database engine (mongodb|sqlite)")
    args = parser.parse_args()

    print "Connecting to database: ", args.db
    if args.dbengine == 'sqlite':
        conn = sqlite3.connect(args.db)
        dbclient = SqliteDB(conn)
    elif args.dbengine == 'mongodb':
        client = pymongo.MongoClient(host=args.host, port=args.port)
        dbclient = MongoDB(client[args.db])
    else:
        raise Exception('please specify a supported database engine !!')
    try:
        finish_corems(dbclient)
    finally:
        dbclient.close()
