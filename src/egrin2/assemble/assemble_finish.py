import argparse
import logging
import sqlite3
import json

import egrin2.query.egrin2_query as e2q
import egrin2.assemble.resample as resample
from egrin2.assemble.assemble_sqlite import SqliteDB


DESCRIPTION = """assemble_finish.py - finish assembling"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None


def __col_resample_pval(dbclient, rows, cols, n_resamples,
                        standardized=True, sig_cutoff=0.05,
                        sort=True, add_override=False, n_jobs=4, keepP=0.1, verbose=False):

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
    to_add = list(filter(lambda col: dbclient.no_col_resamples(col, len(rows), n_resamples), cols))

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

    exp_df = dbclient.find_gene_expressions(rows, cols)
    random_rsd = dbclient.find_col_resamples(len(rows), cols)

    if random_rsd.shape[0] == 0:
        logging.info("Could not find resample DB entry for %d rows in cols %s",
                     len(rows), cols)
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
    pvals.index = pvals.index.values # dbclient.get_col2ids_colout(pvals.index.values)

    if sig_cutoff is not None:
        pvals = pvals[pvals <= sig_cutoff]

    if sort:
        pvals.sort_values()

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
