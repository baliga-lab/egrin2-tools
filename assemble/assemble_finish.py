#!/usr/bin/env

import argparse
import logging
import pymongo
from collections import namedtuple

import query.egrin2_query as e2q
import assemble.resample as resample
import pandas as pd


DESCRIPTION = """assemble_finish.py - finish and dump mongodb"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None


def __check_col_resamples(db, col, n_rows, n_resamples):
    if db.col_resample.find_one({"n_rows": n_rows, "col_id": col, "resamples": {"$gte": n_resamples}}) is None:
        return col


def __col_resample_pval(db, rows, row_type, cols, col_type, n_resamples,
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

    rows_o = rows
    rows = e2q.row2id_batch(db, rows, input_type=row_type)

    if len(rows) == 0:
        logging.info("Please provide an appropriately named array of rows")
        return None

    cols_o = cols
    cols = e2q.col2id_batch(db, cols, input_type=col_type)

    if len(cols) == 0:
        logging.info("Please provide an appropriately named array of cols")
        return None

    # Determine what/how many resamples need to be added to db
    to_add = [__check_col_resamples(db, i, len(rows), n_resamples) for i in cols]
    to_add = [i for i in to_add if i is not None]

    count = 1
    if len(to_add) > 0:
        if add_override:
            logging.info("I need to perform %d random resample(s) of size %d to compute pvals. Please be patient. This may take a while...", len(to_add), n_resamples)
            tmp = resample.col_resample_ind(db, n_rows=len(rows), cols=to_add, n_resamples=n_resamples, keepP=keepP)
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

    exp_df = pd.DataFrame(list(db.gene_expression.find({"col_id": {"$in": cols}, "row_id": {"$in": rows}},
                                                                   { "_id": 0, "col_id": 1, "normalized_expression": 1,
                                                                     "standardized_expression": 1})))
    random_rsd = pd.DataFrame(list(db.col_resample.find({"n_rows": len(rows), "col_id": {"$in": cols}}, {"_id": 0})))

    if random_rsd.shape[0] == 0:
        logging.info("Could not find resample DB entry for %d rows in cols %s", len(rows_o),
                     cols_o)
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
    pvals.index = e2q.col2id_batch(db, pvals.index.values, input_type="col_id", return_field=col_outtype)

    if sig_cutoff is not None:
        pvals = pvals[pvals <= sig_cutoff]

    if sort:
        pvals.sort()

    pvals = pvals.to_frame()
    pvals.columns = ["pval"]

    if pvals.shape[0] == 0:
        logging.info("No cols pass the significance cutoff of %f", sig_cutoff)

    return pvals


def __compute_and_write_col(db, corem, cond_ids, n_resamples=1000):
    logging.info("Adding conditions for corem %d", corem['corem_id'])
    pvals = __col_resample_pval(db, corem['rows'], "row_id", cond_ids, "col_id",
                                n_resamples, keepP=0.05)

    pvals["col_id"] = pvals.index
    d = pvals.to_dict('records')
    print "corem %s -> %s" % (str(corem['_id']), str(d))
    # TODO: enable this line when everything works
    #self.db.corem.update({"_id": corem._id}, {"$set": {"cols": d}})


def finish_corems(db):
    """Finish adding corem info (cols) after resampling. Assumes corem docs already exist"""
    cond_ids = [entry['col_id'] for entry in db["col_info"].find({}, {"_id": 0, "col_id": 1})]

    for corem in db["corem"].find({}, {'_id': 1, "rows": 1, "corem_id": 1}):
        __compute_and_write_col(db, corem, cond_ids)


if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('db', help="MongoDB database name")
    parser.add_argument('--host', default='localhost', required=False, help="MongoDB database host")
    parser.add_argument('--port', default=27017, type=int, required=False, help="MongoDB database port")
    args = parser.parse_args()
    print "Connecting to database: ", args.db
    client = pymongo.MongoClient(host=args.host, port=args.port)
    try:
        finish_corems(client[args.db])
    finally:
        client.close()
