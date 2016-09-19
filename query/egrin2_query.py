#!/usr/bin/env python
"""Tools for querying EGRIN2.0 MongoDB."""

import random
import logging

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import itertools
from bson.code import Code
from collections import defaultdict


def rsd(vals):
    return abs(np.std(vals) / np.mean(vals))


def remove_list_duplicates(l):
    added = set()
    result = []
    for elem in l:
        if elem not in added:
            result.append(elem)
            added.add(elem)
    return result


def find_match(some_row_elem, df, column_name):
    """finds the first row which contains some_row_elem, and returns the
    value specified by the column name"""
    for index, row in df.iterrows():
        if some_row_elem in row.values:
            return row[column_name]
    return None


def row2id_any(db, values, return_field="row_id"):
    """Check name format of rows. If necessary, translate."""
    query = list(db.row_info.find({"$or": [{"row_id": {"$in": values}},
                                           {"egrin2_row_name": {"$in": values}},
                                           {"GI": {"$in": values}},
                                           {"accession": {"$in": values}},
                                           {"name": {"$in": values}},
                                           {"sysName": {"$in": values}}]}))
    if return_field == "all":
        return query
    else:
        return remove_list_duplicates(list([q[return_field] for q in query]))


def row2id_with(db, values, input_field, return_field="row_id"):
    """row2id function searching on a specific input field, which is potentially faster
    than row2id_any()"""
    entry_fields  = {"row_id", "egrin2_row_name", "GI", "accession", "name", "sysName"}
    output_spec = {field: 1 for field in entry_fields}

    if input_field in entry_fields:
        q = {input_field: {"$in": values}}
        query = pd.DataFrame(list(db.row_info.find(q, output_spec)))
        query = query.set_index(input_field)
        to_r = query.loc[values][return_field].tolist()
        return remove_list_duplicates(to_r)
    else:
        raise Exception('row2id_with() called without input_type !!')


def col2id_any(db, values, return_field="col_id"):
    """Check name format of rows. If necessary, translate."""
    query = list(db.col_info.find({"$or": [{"col_id": {"$in": values}},
                                           { "egrin2_col_name": {"$in": values}}]}))
    return [q[return_field] for q in query]


def col2id_with(db, cols, input_field, return_field="col_id"):
    """Check name format of rows. If necessary, translate."""
    if input_field in ["col_id", "egrin2_col_name"]:
        query = pd.DataFrame(list(db.col_info.find({ "$or": [{"col_id": {"$in": cols}},
                                                             {"egrin2_col_name": {"$in": cols}}]},
                                                   {"col_id": 1, "egrin2_col_name": 1})))
        query = query.set_index(input_field)
        return remove_list_duplicates(query.loc[cols][return_field].tolist())
    else:
        raise Exception('col2id_with() called without input type !!!')


# as a first step to clean up the code, move argument type decisions into
# separate functions
def is_row_type(argtype):
    return argtype in {"rows", "row", "gene", "genes"}

def is_col_type(argtype):
     return argtype in {"columns", "column", "col", "cols", "condition", "conditions", "conds"}

def is_gre_type(argtype):
     return argtype in {"motif", "gre", "motc", "motif.gre", "motfs", "gres", "motcs"}

def is_cluster_type(argtype):
     return argtype in {"cluster", "clusters", "bicluster", "biclusters", "bcs"}


def agglom(db, x, x_type, y_type,
           x_input_type=None, y_output_type=None,
           logic="or", gre_lim=10, pval_cutoff=0.05, translate=True):
    """
    Determine enrichment of y given x through bicluster co-membership.
    Available x_type(s)/y_type(s):

    'rows' or 'genes': x should be a gene name or list of genes, eg ["carA","carB"] or [275, 276]
    'columns' or 'conditions': x should be a condition name or list of conditions, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
    'gre':  x should be a GRE ID or list of GRE IDs, eg [4, 19]
    'bicluster': takes or outputs bicluster '_id'
    """
    logging.info('Using "%s" logic', logic)

    def compute_p(i, M, N):
        z = i.counts # n black balls in draw
        n = i.all_counts # num black balls tot
        M = M
        # M = M - n  # n white balls
        N = N # num drawn
        prb =  hypergeom.sf(z, M, n, N)
        return prb

    if x_type is None:
        logging.info("""Please supply an x_type for your query.
Types include: 'rows' (genes), 'columns' (conditions), 'gres'""")
    if y_type is None:
        logging.info("""Please supply an y_type for your query.
Types include: 'rows' (genes), 'columns' (conditions), 'gres'. Biclusters will be returned by default.""")
        y_type = "cluster"

    if type(x) == str or type(x) == int:
        x = [x]  # single

    # Check input types
    if is_row_type(x_type):
        x_type = "rows"
        x_o = x
        if x_input_type is not None:
            x = row2id_with(db, x, x_input_type, return_field="row_id")
        else:
            x = row2id_any(db, x, return_field="row_id")

        if len(x) == 0:
            logging.info("Cannot translate row names: %s", x_o)
            return None

    elif is_col_type(x_type):
        x_type = "columns"
        x_o = x
        if x_input_type is not None:
            x = col2id_with(db, x, x_input_type, return_field="col_id")
        else:
            x = col2id_any(db, x, return_field="col_id")

        if len(x) == 0:
            logging.info("Cannot translate col names: %s", x_o)
            return None

    elif is_gre_type(x_type):
        x_type = "gre_id"

    elif is_cluster_type(x_type):
        logging.warn("I hope you are using cluster '_id'!!! Otherwise the results might surprise you...")
        x_type = "_id"

    else:
        logging.error("Can't recognize your 'x_type' argument.")
        return None

    # Check output types
    if is_row_type(y_type):
        y_type = "rows"
    elif is_col_type(y_type):
        y_type = "columns"
    elif is_gre_type(y_type):
        y_type = "gre_id"
    elif is_cluster_type(y_type):
        logging.warn("Will return bicluster _id. The results might surprise you...")
        y_type = "_id"
    else:
        logging.error("Can't recognize your 'y_type' argument.")
        return None

    # Compose query
    if logic in ["and","or","nor"]:
        q = {"$" + logic: [{x_type: i} for i in x]}
        o = {y_type: 1}

        if x_type == "gre_id":
            queryPre = pd.DataFrame(list(db.motif_info.find(q, {"cluster_id" : 1})))["cluster_id"].tolist()
            query = pd.DataFrame(list(db.bicluster_info.find({"_id": {"$in": queryPre}}, o)))

        elif y_type == "gre_id":
            queryPre = pd.DataFrame(list(db.bicluster_info.find(q, {"_id" : 1})))["_id"].tolist()
            query = pd.DataFrame(list(db.motif_info.find({"cluster_id": {"$in": queryPre}}, {y_type: 1})))

        else:
            query = pd.DataFrame(list(db.bicluster_info.find(q, o)))
    else:
        logging.error("I don't recognize the logic you are trying to use. 'logic' must be 'and', 'or', or 'nor'.")
        return None

    if query.shape[0] > 0:

        mapColumns = Code("""
            function () {
            this.columns.forEach(function(z) {
                emit(z, 1);
                });
            }
            """)
        mapRows = Code("""
            function () {
            this.rows.forEach(function(z) {
                emit(z, 1);
                });
            }
            """)
        mapGREs = Code("""
            function () {
            print(this);
            emit(this.gre_id, 1);
            }
            """)
        reduce = Code("""
            function (key, values) {
                     var total = 0;
                     for (var i = 0; i < values.length; i++) {
                    total += values[i];
                    }
                return total;
                 }
                 """)

        if y_type == "_id":
            return query
        else:
            if y_type == "rows":

                if db.rowsCount_mapreduce.count() == 0:
                    logging.info("Initializing MapReduce lookup table. Future queries will be much faster!")
                    db.bicluster_info.map_reduce(mapRows, reduce, "rowsCount_mapreduce")
                else:
                    # do spot check to make sure mapreduce is up to date
                    random_id = random.randint(0, db.rowsCount_mapreduce.count())
                    ref = db.rowsCount_mapreduce.find_one({"_id": random_id})["value"]
                    test = db.bicluster_info.find({"rows": random_id}).count()

                    if ref != test:
                        logging.info("Initializing MapReduce lookup table. Future queries will be much faster!")
                        db.bicluster_info.map_reduce(mapRows, reduce, "rowsCount_mapreduce")

                rows = pd.Series(list(itertools.chain(*query.rows.tolist()))).value_counts().to_frame("counts")

                # filter out rows that aren't in the database - i.e. not annotated in MicrobesOnline
                in_db = pd.DataFrame(list(db.row_info.find({}, {"_id": 0, "row_id": 1}))).row_id.tolist()
                common_rows = list(set(rows.index).intersection(set(in_db)))
                rows = rows.loc[common_rows]

                # find all bicluster counts
                all_counts = pd.DataFrame(list(db.rowsCount_mapreduce.find()))
                all_counts = all_counts.set_index("_id")

                # combine two data frames
                to_r = rows.join(all_counts).sort_values("counts", ascending=False)
                to_r.columns = ["counts","all_counts"]
                to_r = to_r.drop_duplicates()

                if translate:
                    to_r.index = row2id_with(db, to_r.index.tolist(), "row_id", return_field="egrin2_row_name")

            if y_type == "columns":

                if db.columnsCount_mapreduce.count() == 0:
                    logging.info("Initializing MapReduce lookup table. Future queries will be much faster!")
                    db.bicluster_info.map_reduce(mapColumns, reduce, "columnsCount_mapreduce")

                else:
                    # do spot check to make sure mapreduce is up to date
                    random_id = random.randint(0, db.columnsCount_mapreduce.count())
                    ref = db.columnsCount_mapreduce.find_one({"_id": random_id})["value"]
                    test = db.bicluster_info.find({"columns": random_id}).count()

                    if ref != test:
                        logging.info("Initializing MapReduce lookup table. Future queries will be much faster!")
                        db.bicluster_info.map_reduce(mapColumns, reduce, "columnsCount_mapreduce")

                if db.columnsCount_mapreduce.count() == 0:
                    db.bicluster_info.map_reduce(mapColumns, reduce, "columnsCount_mapreduce")

                cols = pd.Series(list(itertools.chain(*query["columns"].tolist()))).value_counts().to_frame("counts")

                # filter out columns that aren't in the database - i.e. not annotated in MicrobesOnline
                in_db = pd.DataFrame(list(db.col_info.find({}, {"_id": 0, "col_id": 1}))).col_id.tolist()
                common_cols = list(set(cols.index).intersection(set(in_db)))
                cols = cols.loc[common_cols]

                # find all bicluster counts
                all_counts = pd.DataFrame(list(db.columnsCount_mapreduce.find()))
                all_counts = all_counts.set_index("_id")

                # combine two data frames
                to_r = cols.join(all_counts).sort_values("counts", ascending=False)
                to_r.columns = ["counts","all_counts"]

                if translate:
                    to_r.index = col2id_with(db, to_r.index.tolist(), "col_id", return_field="egrin2_col_name")

            if y_type == "gre_id":

                if db.gresCount_mapreduce.count() == 0:
                    logging.info("Initializing MapReduce lookup table. Future queries will be much faster!")
                    db.motif_info.map_reduce(mapGREs, reduce, "gresCount_mapreduce")

                else:
                    # do spot check to make sure mapreduce is up to date
                    random_id = random.randint(0, db.gresCount_mapreduce.count())
                    ref = db.gresCount_mapreduce.find_one({"_id": random_id})["value"]
                    test = db.motif_info.find({"gre_id": random_id}).count()
                    if ref != test:
                        logging.info("Initializing MapReduce lookup table. Future queries will be much faster!")
                        db.motif_info.map_reduce(mapGREs,reduce,"gresCount_mapreduce")

                gres = query.gre_id.tolist()
                gres = list(filter(lambda x: x != "NaN", gres))
                gres = pd.Series(gres).value_counts().to_frame("counts")

                # find all bicluster counts
                all_counts = pd.DataFrame(list(db.gresCount_mapreduce.find()))
                all_counts = all_counts.set_index("_id")

                # combine two data frames
                to_r = gres.join(all_counts).sort_values("counts", ascending=False)
                to_r.columns = ["counts","all_counts"]

                # filter by GREs with more than 10 instances
                to_r = to_r.loc[to_r.all_counts>=gre_lim, :]

            to_r["pval"] = to_r.apply(compute_p, axis=1, M=db.bicluster_info.count(), N=query.shape[0])
            to_r["qval_BH"] = multipletests(to_r.pval, method='fdr_bh')[1]
            to_r["qval_bonferroni"] = multipletests(to_r.pval, method='bonferroni')[1]
            to_r = to_r.sort_values(["pval","counts"], ascending=True)

            # only return below pval cutoff
            to_r = to_r.loc[to_r.pval <= pval_cutoff, :]
            return to_r

    else:
        logging.info("Could not find any biclusters matching your criteria")
        return None


def find_gre_positions(db, gre_ids):
    """return a list of all positions of motifs associated with a list of GREs
    return a dataframe of (gre_id, start, stop) rows
    """
    positions = []
    cluster_mot_pos = {}

    gres =  db.motif_info.find({'gre_id': { '$in': gre_ids}}, {'_id': 0, 'gre_id': 1, 'cluster_id': 1, 'motif_num': 1})
    for entry in gres:
        gre_id = entry['gre_id']
        cluster_id = entry['cluster_id']
        motif_num = entry['motif_num']
        key = '%s:%d' % (str(cluster_id), motif_num)
        if key in cluster_mot_pos:
            mot_pos = cluster_mot_pos[key]
        else:
            mot_pos = []
            for entry in db.fimo.find({'cluster_id': cluster_id, 'motif_num': motif_num}, {'start': 1, 'stop': 1}):
                mot_pos.append((entry['start'], entry['stop']))
        for start, stop in mot_pos:
            positions.append((gre_id, start, stop))
    return pd.DataFrame(positions, columns=['gre_id', 'start', 'stop'])


def find_fimo(db, start=None, stop=None, locusId=None, strand=None, mot_pval_cutoff=None, filterby=None, filter_type=None,
              filterby_input_type=None, use_fimo_small=True,
              logic="or", return_format="file", outfile=None, tosingle=True):
    """
    Find motifs/GREs that fall within a specific range. Filter by biclusters/genes/conditions/etc. Optionally write to GGBweb compatible .gsif format.
    IMPORTANT!!!! Only filtering by GREs currently supported

    Parameters:

    -- db: mongo database object
    -- start: genomic start, chromosome start if None
    -- stop: genomic start, chromosome end if None
    -- locusId: chromosome Id, (i.e. MicrobesOnline scaffoldId or NCBI_RefSeq)
    -- strand: fetch stand specific predictions (not supported currently)
    -- mot_pval_cutoff: only retrieve fimo scan matches below this pval_cutoff
    -- filterby: elements through which to filter fimo scans, e.g. gres
    -- filter_type: type of filterby, e.g. gres
    -- filterby_input_type: name format of filterby. Speeds up name translation
    -- use_fimo_small: use fimo_small collection (highly significant matches only). Speed up query.
    -- logic: logcal operation to apply to filterby, ie and, or, nor
    -- return_format: how should the fimo tracks be returned, only "file" currently support (for upload to GGBweb)
    -- outfile: path and name of output file in .gsif format
    -- tosingle: return one filter per filterby object (eg GRE) or all in a single file (for early GGBweb dev support)
    """
    def getBCs(x, x_type):
        if x is None:
            to_r = pd.DataFrame(list(db.motif_info.find({}, {"cluster_id": 1, "gre_id": 1})))
        else:
            to_r = pd.DataFrame(list(db.motif_info.find({x_type: x}, {"cluster_id": 1, "gre_id": 1})))
        return(to_r.loc[:, ["gre_id" ,"cluster_id"]])

    def aggSeq(x):
        def count(y):
            return range(y.start, y.stop + 1)

        if x.shape[0] == 0:
            return []

        else:
            to_r = [count(x.iloc[i]) for i in range(x.shape[0])]
            to_r = list(itertools.chain(*to_r))
            to_r.sort_values()

        return to_r

    db_chr = pd.DataFrame(list(db.genome.find({}, {"scaffoldId": 1, "NCBI_RefSeq": 1})))
    db_scaffoldId = db_chr.scaffoldId.tolist()
    db_NCBI_RefSeq = db_chr.NCBI_RefSeq.tolist()

    if locusId is None:
        logging.error("""Please provide a chromosome Locus ID.
This is probably the scaffoldID from MicrobesOnline.
LocusIds in database %s include:

  - ScaffoldId
  - %s
  - NCBI_RefSeq
  - %s""",
                      db.name, (", ").join(db_scaffoldId), (", ").join(db_NCBI_RefSeq))
        return None

    locusId = str(locusId)

    if locusId not in db_scaffoldId and locusId not in db_NCBI_RefSeq:
        logging.error("""LocusId %s not in EGRIN 2.0 database %s.
LocusIds in this database include:

  - ScaffoldId
  - %s
  - NCBI_RefSeq
  - %s""",
                      locusId, db.name,
                      (", ").join(db_scaffoldId),
                      (", ").join(db_NCBI_RefSeq))
        return None

    chromosome = db.genome.find_one({"$or": [{"scaffoldId": locusId}, {"NCBI_RefSeq": locusId}]})
    scaffoldId = chromosome["scaffoldId"]
    ncbi = chromosome["NCBI_RefSeq"]

    # if start/stop is None, assume whole chromosome
    if start is None:
        logging.info("Start not provided. Assuming beginning of chromosome")
        start = 0

    if stop is None:
        logging.info("Stop not provided. Assuming end of chromosome")
        stop = len( chromosome["sequence"] )

    if use_fimo_small:
        fimo_collection = "fimo_small"
    else:
        fimo_collection = "fimo"

    if filterby is None:
        logging.info("No filter applied")

    if filter_type is not None:
        logging.warn("Many of these filters are not supported currently. Only GREs!!!")
        if is_row_type(filter_type):
            filter_type = "rows"
            filterby_o = filterby
            filterby = row2id_with(db, filterby, filter_input_type, return_field="row_id")
            filterby = list(set(filterby))
            if len(filterby) == 0:
                logging.error("Cannot translate row names: %s", filterby_o)
                return None

        elif is_col_type(filter_type):
            filter_type = "columns"
            filterby_o = filterby
            filterby = col2id_with(db, filterby, filterby_input_type, return_field="col_id")
            filterby = list(set(filterby))

            if len(filterby) == 0:
                logging.error("Cannot translate col names: %s", filterby_o)
                return None

        elif is_gre_type(filter_type):
            filter_type = "gre_id"

        elif is_cluster_type(filter_type):
            logging.warn("I hope you are using cluster '_id'!!! Otherwise the results might surprise you...")
            filter_type = "_id"

        logging.info("Filtering motifs by %s", filter_type)
        bcs_df= pd.concat([getBCs(i, filter_type) for i in filterby], ignore_index=True)
        mots = pd.DataFrame(list(db[fimo_collection].find({"start": {"$gte": start},
                                                           "stop": {"$lte": stop },
                                                           "cluster_id": {"$in": bcs_df.cluster_id.tolist()},
                                                           "scaffoldId": scaffoldId})))
        mots = pd.merge(mots, bcs_df, on="cluster_id")
    else:
        bcs_df= pd.concat([getBCs(None, filter_type)], ignore_index = True )
        mots = pd.DataFrame(list(db[fimo_collection].find({ "start": {"$gte": start}, "stop": {"$lte": stop}, "scaffoldId": scaffoldId})))
        mots = pd.merge(mots, bcs_df, on="cluster_id")

    gre_scans = mots.groupby("gre_id").apply(aggSeq)

    if return_format == "dictionary":
        gre_scans = {i: pd.Series( gre_scans.loc[i]).value_counts().sort_index() for i in gre_scans.index}

    if return_format == "file":
        # return with file format ready to save for GGBweb
        # start stop strand chr value id
        tmp_gre_scans = []

        for i in gre_scans.index:
            tmp_df = pd.DataFrame(columns=["start","end","strand","chr","value","id"])
            gre_counts = pd.Series( gre_scans.loc[i] ).value_counts().sort_index()
            tmp_df["start"] =  gre_counts.index
            tmp_df["end"] =  gre_counts.index
            tmp_df["strand"] = "+"
            tmp_df["chr"] = ncbi
            tmp_df["value"] = gre_counts.values
            tmp_df["id"] = "GRE_" + str(i)
            tmp_gre_scans.append(tmp_df)

        gre_scans = pd.concat(tmp_gre_scans, ignore_index=True)

        if outfile is not None:
            def dfsave(df, fname):
                fname = list(set(df.id))[0] + "_" + fname
                logging.info("Writing file '%s'", fname)

                df = df.drop('id', 1)
                df.to_csv(fname, sep="\t", index=False)
                return None

            if tosingle:
                gre_scans_grouped = gre_scans.groupby("id")
                tmp = [dfsave(gre_scans_grouped.get_group(i), fname=outfile)
                       for i in gre_scans_grouped.groups.keys()]
                return None
            else:
                gres_scans.to_csv(outfile, sep="\t", index=False)
                return None

    return gre_scans


def find_corem_info(db, x, x_type="corem_id", x_input_type=None, y_type="genes", y_return_field=None,
                    count=False, logic="or"):

    """
    Fetch corem-related info 'y' given query 'x'.

    Available x_type(s)/y_type(s):

    'corem_id': x is corem_id (int), eg [1]
    'rows' or 'genes': x should be a gene name or list of genes, eg ["carA","carB"] or [275, 276]
    'columns' or 'conditions': x should be a condition or list of conditions, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
    'gre': x should be a GRE or a list of GRE IDs, eg [4, 19]
    'edge': x should be an edge or list of edges. Genes in edges should be separated by '-', eg ["carA-carB"] or ["275-276"]
    """
    TYPE_MAP = {'rows': 'rows', 'row': 'rows', 'gene': 'rows', 'genes': 'rows',
                'columns': 'cols.col_id', 'column': 'cols.col_id', 'col': 'cols.col_id', 'cols': 'cols.col_id',
                'condition': 'cols.col_id', 'conditions': 'cols.col_id', 'conds': 'cols.col_id',
                'motif': 'gre_id', 'gre': 'gre_id', 'motc': 'gre_id', 'motif.gre': 'gre_id', 'motifs': 'gre_id',
                'gres': 'gre_id', 'motcs': 'gre_id',
                'corem_id': 'corem_id', 'corem': 'corem_id', 'corems': 'corem_id',
                'edge': 'edges', 'edges': 'edges'}

    if x_type is None:
        logging.info("Please supply an x_type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres', edges")
    if y_type is None:
        logging.info("Please supply an y_type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres', edges")

    if type(x) == str or type(x) == int:
        x = [x]  # single

    # Check input types
    x_type = TYPE_MAP[x_type]
    if x_type == "rows":
        x_o = x
        x = row2id_with(db, x, x_input_type, return_field="row_id")
        x = list(set(x))
        if len(x) == 0:
            logging.error("Cannot translate row names: %s", x_o)
            return None
    elif x_type == "cols.col_id":
        x_o = x
        x = col2id_with(db, x, x_input_type, return_field="col_id")
        x = list(set(x))
        if len(x) == 0:
            logging.error("Cannot translate row names: %s", x_o)
            return None
    elif x_type == "edges":
        x_new = []
        for i in x:
            i_trans = row2id_with(db, i.split("-"), x_input_type, return_field="row_id")
            i_trans = [str(j) for j in i_trans]
            x_new.append("-".join(i_trans))
            i_trans.reverse()
            x_new.append("-".join(i_trans))

        x = x_new

        if len(x) == 0:
            logging.error("Cannot translate row names: %s", x_o)
            return None

    y_type_original = y_type
    y_type = TYPE_MAP[y_type]

    if logic in {"and","or","nor"}:
        if logic == "and" and x_type == "corem_id":
            q = {"$or": [{x_type: i} for i in x]}

        else:
            q = {"$" + logic: [{x_type: i} for i in x]}

        o = {y_type: 1}
        query = pd.DataFrame(list(db.corem.find(q, o)))

    else:
        logging.error("I don't recognize the logic you are trying to use. 'logic' must be 'and', 'or', or 'nor'.")
        return None

    if y_type == "rows":
        if y_return_field is None:
            y_return_field = "egrin2_row_name"

        if query.shape[0] > 1:
            to_r = list(itertools.chain( *query.rows.values.tolist()))

            if logic == "and":
                to_r = pd.Series(to_r).value_counts()

                if count:
                    to_r = to_r[to_r >= query.shape[0]]

                    if to_r.shape[0] > 0:
                        to_r.index = row2id_with(db, to_r.index.tolist(), "row_id", return_field=y_return_field)
                    else:
                        logging.error("No genes found")
                        return None
                else:
                    to_r = to_r[to_r >= query.shape[0]].index.tolist()

                    if len(to_r):
                        to_r = row2id_with(db, to_r, "row_id", return_field=y_return_field)
                        to_r.sort_values()
                    else:
                        logging.error("No genes found")
                        return None
            else:
                if count:
                    to_r = pd.Series(to_r).value_counts()

                    if to_r.shape[0] > 0:
                        to_r.index = row2id_with(db, to_r.index.tolist(), "row_id", return_field=y_return_field)
                    else:
                        logging.error("No genes found")
                        return None
                else:
                    to_r = list(set(to_r))
                    if len(to_r):
                        to_r = row2id_with(db, to_r, "row_id", return_field=y_return_field)
                        to_r.sort_values()
                    else:
                        logging.error("No genes found")
                        return None
        else:
            to_r = row2id_with(db, query.rows[0], "row_id", return_field=y_return_field)

    elif y_type == "cols.col_id":
        if y_return_field is None:
            y_return_field = "egrin2_col_name"

        if query.shape[0] > 1:
            to_r = [int(i["col_id"]) for i in list(itertools.chain(*query.cols.values.tolist())) if i["col_id"] if type(i["col_id"]) is float]

            if logic == "and":
                to_r = pd.Series(to_r).value_counts()

                if count:
                    to_r = to_r[to_r >= query.shape[0]]
                    if to_r.shape[0] > 0:
                        to_r.index = col2id_with(db, to_r.index.tolist(), "col_id", return_field=y_return_field)
                    else:
                        logging.error("No conditions found")
                        return None
                else:
                    to_r = to_r[to_r >= query.shape[0]].index.tolist()

                    if len(to_r):
                        to_r = col2id_with(db, to_r, 'col_id', return_field=y_return_field)
                        to_r.sort_values()
                    else:
                        logging.error("No conditions found")
                        return None
            else:
                if count:
                    to_r = pd.Series(to_r).value_counts()
                    if to_r.shape[0] > 0:
                        to_r.index = col2id_with(db, to_r.index.tolist(), 'col_id', return_field=y_return_field)
                    else:
                        logging.error("No conditions found")
                        return None
                else:
                    to_r = list(set(to_r))
                    if len(to_r):
                        to_r = col2id_with(db, to_r, 'col_id', return_field=y_return_field)
                        to_r.sort_values()
                    else:
                        logging.error("No conditions found")
                        return None
        else:
            to_r = [int(i["col_id"]) for i in list(itertools.chain( *query.cols.values.tolist())) if i["col_id"] if type(i["col_id"]) is float]
            to_r = col2id_with(db, to_r, 'col_id', return_field=y_return_field)

    elif y_type == "corem_id":
        if query.shape[0] > 0:
            to_r = query.corem_id.tolist()
        else:
            to_r = None

    elif y_type == "edges":
        if y_return_field is None:
            y_return_field = "egrin2_row_name"
        to_r = list(itertools.chain(*query.edges.values.tolist()))
        to_r_new = []

        for i in to_r:
            i_trans = row2id_with(db, [int(j) for j in i.split("-")], "row_id", return_field=y_return_field)

            i_trans = [str(j) for j in i_trans]
            to_r_new.append("-".join(i_trans))
        to_r = to_r_new
        to_r.sort_values()
    elif y_type == "gre_id":
        """TODO: GRE's for a corem"""
        logging.warn("GREs detection is not currently supported for corems. Please use the `agglom` function to find GREs enriched in biclusters containing corem genes instead.")
        to_r = None
    else:
        logging.error("Could not find corems matching your query")
        to_r = None

    if to_r is not None:
        to_r = pd.DataFrame(to_r)
        to_r.columns = [y_type_original]

    return to_r


def find_gene_expression(db, rows=None, cols=None, standardized=True):
    """
    Fetch gene expression given rows and columns.

    Parameters:
    -- db: mongo database object
    -- rows: list of rows/genes in format recognized by row2id_with() (i.e. some name present in MicrobesOnline)
    -- cols: list of columns/conditions in format recognized by col2id_with (i.e. name in ratios matrix)
    -- standardized: fetch standardized data if True, otherwise raw (normalized) data

    """
    input_type_rows = None
    input_type_cols = None

    if rows is None:
        # assume all genes
        rows = pd.DataFrame(list(db.row_info.find({}, {"row_id":1}))).row_id.tolist()
        input_type_rows = "row_id"

    if cols is None:
        # assume all cols
        cols = pd.DataFrame(list(db.col_info.find({}, {"col_id":1}))).col_id.tolist()
        input_type_cols = "col_id"

    if type(rows) == str or type(rows) == int:
        rows = [rows]  # single

    if type(cols) == str or type(cols) == int:
        cols = [cols]  # single

    # translate rows/cols
    rows = row2id_with(db, rows, input_type_rows, return_field="row_id")
    cols = col2id_with(db, cols, input_type_cols, return_field="col_id")

    if len(rows) > 1000 or len(cols) > 1000:
        logging.warn("This is a large query. Please be patient. If you need faster access, I would suggest saving this matrix and loading directly from file.")

    # get expression data
    data = pd.DataFrame(None, columns=cols, index=rows)
    query = pd.DataFrame (list(db.gene_expression.find({"$and": [{"row_id": {"$in": rows}}, {"col_id": {"$in": cols}}]})))

    for i in query:
        if standardized:
            data = query.pivot(index="row_id", columns="col_id", values="standardized_expression")
        else:
            data = query.pivot(index="row_id", columns="col_id", values="raw_expression")

    data.index = row2id_with(db, data.index.tolist(), "row_id", return_field="egrin2_row_name")
    data.columns = col2id_with(db, data.columns.tolist(), 'col_id', return_field="egrin2_col_name")
    data = data.sort_index()
    data = data.reindex_axis(sorted(data.columns), axis=1)

    return data


def export_ggbweb(db, genes=None, outfile=None):
    """
    Write gene module in GGBweb .gsif format

    Parameters:
    -- db: mongo database object
    -- genes: list of genes in format recognized by row2id_with() (some name present in MicrobesOnline)
    -- outfile: path/name of module file to write
    """
    if genes is None:
        logging.error("Please provide a gene or list of genes")
        return None

    if type(genes) == str or type(genes) == int:
        genes = [genes] # single

    def locFormat(x, gene):
        return pd.DataFrame([gene, "E.colik-12", str(x.start), str(x.stop)])

    gene_info = pd.DataFrame(row2id_any(db, genes, return_field="all"))

    to_r = pd.concat([locFormat(gene_info.iloc[i], gene_info.loc[i, "egrin2_row_name"]) for i in range(gene_info.shape[0])], 2).T

    if outfile is not None:
        logging.info("Module written to: '%s'", os.path.abspath(outfile))
        to_r.to_csv(os.path.abspath( outfile ), sep="\t", index=False, header=False)
        return None

    return to_r


def find_motifs(db, x, x_type, output_type=["data_frame", "array"][0]):
    """
    Find, retrieve, and format motif PWM

    Parameters:
    -- db: mongo database object
    -- x: list of motifs to retrieve. Can be combined name separated by `_`, eg to return motif #1 from bicluster #10 in run #3 supply x = "3_10_1". See below for x_type.
    -- x_type: format of x, e.g. gre. Can be combined name separated by `_`, eg to return motif #1 from bicluster #10 in run #3 supply x_type = "run_cluster_motif"
    -- output_type: format of output. `data_frame` for pretty printing, `array` for viz by weblogo

    Returns DataFrame (or list of DataFrames) containing PWM (more precisely PPMs)
    """
    def reshapePWM(x):
        """
        x is PWM list from MongoDB to be converted to DataFrame
        """
        to_r = pd.DataFrame(x)
        to_r = to_r.set_index("row").sort_index()

        # make sure order is consistent with weblogo
        to_r = to_r.loc[:, ["a","c","g","t"]]
        return to_r

    if type(x) == str or type(x) == int:
        x = [x] # single

    if len(x_type.split("_")) > 1:

        x_type_split = x_type.split("_")
        logging.info("You have indicated a combined name type: %s.\n\nIf this is not your intent, I suggest you remove the `_`",
                     (", ").join(x_type_split))

        # change input names
        for i in range(len(x_type_split)):
            if x_type_split[i] == "run" or x_type_split[i] == "runs" :
                x_type_split[i] = "run_id"
            elif x_type_split[i] == "bicluster" or x_type_split[i] == "cluster" :
                x_type_split[i] = "cluster"
            elif x_type_split[i] == "motif" or x_type_split[i] == "gre" :
                x_type_split[i] = "motif_num"
            else:
                logging.error("Could not match x_type: %s", x_type_split[i])

        # format query
        q1 = {}
        q1[ "$or" ] = []
        q2 = pd.DataFrame()

        # makes sure lengths of x and x_type are the same
        for i in x:
            i_split = [int(j) for j in i.split("_")]

            if len(i_split) != len(x_type_split):
                logging.error("Number of combined x_types does not match x")
            else:
                i_dict = dict(zip(x_type_split, i_split))
                q2 = pd.concat([q2, pd.DataFrame.from_dict(i_dict,orient="index" ).T])
                if "motif_num" in x_type_split:
                    i_dict.pop("motif_num")

                q1["$or"].append(i_dict)

        # doesn't matter that "motif_info" is here
        o = {i: 1 for i in x_type_split}

        # get biclusters
        prequery = pd.DataFrame(list(db.bicluster_info.find(q1, o)))
        prequery = prequery.merge(q2, on=list(set(prequery.columns.tolist()) & set(x_type_split)))
        prequery = prequery.rename(columns={'_id' : 'cluster_id'})

        if "motif_num" in x_type_split:
            # specific bc motifs requested
            q = {"$or": prequery.loc[:, ["cluster_id","motif_num"]].to_dict("records")}
        else:
            q = {"$or": prequery.loc[:, ["cluster_id"]].to_dict("records")}

        query = pd.DataFrame(list(db.motif_info.find(q)))
        query = query.merge(prequery, on=list(set(query.columns.tolist()) & set(prequery.columns.tolist())))

        query["x_type"] = ""
        for i in range(query.shape[0]):
            query.loc[i, "x_type"] = ("_").join([str(j) for j in query.loc[i, x_type_split].tolist()])

        query = query.set_index("x_type")
        # REORDER!!!
        query = query.loc[x, :]

    elif is_gre_type(x_type):
        x_type = "gre_id"
        q = {x_type: {"$in": x}}
        o = {"_id": 0, "pwm": 1}
        query = pd.DataFrame(list(db.corem.find(q, o)))
    else:
        logging.error("Cannot recognize x_type = %s", x_type)

    if output_type == "array":
        to_r = {query.index[i]: reshapePWM(query.iloc[i].pwm).as_matrix() for i in range(query.shape[0])}
    else:
        to_r = {query.index[i]: reshapePWM(query.iloc[i].pwm) for i in range( query.shape[0])}

    if len(to_r) == 1:
        to_r = to_r.values()[0]

    return to_r


def gre_motifs(db, gre_id, evalue=None):
    """Returns a list of PSSMs for the given GRE id.
    If evalue is specified, only the motifs with a smaller evalue are returned
    """
    result = []
    query = {'gre_id': gre_id}
    if evalue is not None:
        query['evalue'] = {'$lt': evalue}

    pwms = db.motif_info.find(query, {'pwm': 1, '_id': 0})

    for i, pwm in enumerate(pwms):
        rows = pwm['pwm']
        result.append([[row['a'], row['g'], row['c'], row['t']] for row in rows])
    return result


def gre_motifs_batch(db, gre_ids, evalue=None):
    """Returns a list of PSSMs for the given GRE id.
    If evalue is specified, only the motifs with a smaller evalue are returned
    """
    result = defaultdict(list)
    query = {'gre_id': {"$in": gre_ids}}
    if evalue is not None:
        query['evalue'] = {'$lt': evalue}

    gre_pwms = db.motif_info.find(query, {'gre_id': 1, 'pwm': 1, '_id': 0})

    for gre_pwm in gre_pwms:
        result[gre_pwm['gre_id']].append([[row['a'], row['g'], row['c'], row['t']] for row in gre_pwm['pwm']])
    return result
