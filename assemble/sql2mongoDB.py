#!/usr/bin/env python

"""
Initialize  mongo database from individual cMonkey
runs (SQLite) plus some additional tables and run_id column,
including motifs and motif clusters
"""
import os
import datetime
import glob
import sys
import gzip
import time
from urllib2 import urlopen, URLError, HTTPError
from zipfile import ZipFile
import itertools
import logging

import numpy as np
import pandas as pd

import sqlite3
import pymongo
from bson.objectid import ObjectId
from Bio import SeqIO


DEFAULT_MOTIF_OUT_FILENAME = "out.mot_metaclustering.txt.I45.txt"


###########################################################################
### Helper functions
###
### These module private functions support the database
### import
###########################################################################


def _gre2motif_path(ensembledir, gre2motif):
    """Determine GRE clustering file path"""
    if gre2motif is None:
        motif_out_path = os.path.join(ensembledir, DEFAULT_MOTIF_OUT_FILENAME)

        if os.path.isfile(motif_out_path):
            return motif_out_path
        else:
            logging.warn("I cannot find a GRE clustering file. If you want to assign GREs, please specify this file.")
            return None
    else:
        return gre2motif


def _get_ncbi_code(ncbi_code, resultdb_path):
    """Determine the NCBI code either from the ncbi_code argument or the given result database"""
    if ncbi_code is None:
        # try to find ncbi_code in cmonkey_run.db
        conn = sqlite3.connect(resultdb_path)
        c = conn.cursor()
        try:
            c.execute("SELECT ncbi_code FROM run_infos")
            return c.fetchone()[0]
        except sqlite3.Error as e:
            logging.warn("Could not find NCBI Genome ID in cmonkey_run.db: %s", resultdb_path)
        finally:
            c.close()
            conn.close()
    else:
        return int(ncbi_code)


def _available_cmonkey_resultdb_files(ensembledir, prefix):
    """Search the ensemble directory for available cmonkey runs and return
    a list of the paths of cmonkey runs found"""
    result = sorted(glob.glob(os.path.join(ensembledir, "%s???/cmonkey_run.db" % prefix)))
    if len(result) ==  0:
        logging.warn("""I cannot find any cMonkey SQLite databases in the current directory: %s
Make sure 'ensembledir' variable points to the location of your cMonkey-2 ensemble results.""",
                     os.getcwd())
    return result


def _download_url(url, save_name='tmp', num_retries=5):
    """Download the file at the specified location"""
    count = 1
    while count < num_retries:
        try:
            f = urlopen(url)
            logging.info("downloading '%s'", url)

            # Open our local file for writing
            with open(os.path.basename(save_name), "wb") as local_file:
              local_file.write(f.read())
            break
        # handle errors
        except HTTPError, e:
            logging.error("HTTP Error, code: %s URL: %s", str(e.code), url)
            logging.info("Trying to connect again. Attempt %d of 5", count)
            count += 1
        except URLError, e:
            logging.error("URL Error, Reason: %s URL: %s", e.reason, url)
            logging.info("Trying to connect again. Attempt %d of 5", count)
            count += 1
    return count < num_retries


def _valid_cmonkey_results(db, resultdb_files, db_run_override):
    """make sure the runs have data!!!"""
    def keep_this_run(dbpath):
        """if db_run_override is false, only include runs that are not yet in the database"""
        run_name = dbpath.split("/")[-2]
        return db_run_override or db.ensemble_info.find({"run_name": run_name}).count() == 0

    result = []
    ensemble_info_collection = db.ensemble_info
    for dbpath in resultdb_files:
        conn = sqlite3.connect(dbpath)
        c = conn.cursor()
        try:
            c.execute("select finish_time from run_infos")
            row = c.fetchone()
            if row is not None and row[0] is not None and keep_this_run(dbpath):
                result.append(dbpath)
            else:
                if row is None or row[0] is None:
                    logging.warn("incomplete run: '%s'", dbpath)
                else:
                    logging.warn("already in the database: '%s'", dbpath)
        except Exception:
            logging.warn("incomplete run: '%s'", dbpath)
        finally:
            c.close()
            conn.close()
    return result


def _not_exists(collection, doc):
    """check if the document exists in the specified MongoDB collection"""
    return collection.find(doc).count() == 0


def _not_exists_key_value(collection, key, value):
    """check if the document exists in the specified MongoDB collection
    using a key-value pair"""
    return collection.find({key: value}).count() == 0

def _not_exists_key_value2(collection, key1, value1, key2, value2):
    """check if the document exists in the specified MongoDB collection
    using two key-value pairs"""
    return collection.find({key1: value1, key2: value2}).count() == 0


def _get_run2id(db, db_files):
    """make run2id"""
    run_names_and_ids = [(run_info["run_name"], run_info["run_id"])
                         for run_info in db.ensemble_info.find()]
    run_names = [rni[0] for rni in run_names_and_ids]
    run_ids = [rni[1] for rni in run_names_and_ids]

    for dbpath in db_files:
        run_name = dbpath.split("/")[-2]
        if run_name not in run_names:
            run_names.append(run_name)
            if len(run_ids) > 0:
                run_ids.append(max(run_ids) + 1)
            else:
                run_ids.append(0)

    return pd.DataFrame(zip(run_ids, run_names), index=run_names, columns=["run_id", "run_name"])


def _import_genome(db, genome_file, ncbi_code):
    """Import genome information into database"""
    logging.info("Downloading genome information for NCBI taxonomy ID: %s", ncbi_code)
    if genome_file is None:
        logging.info("No custom genome annotation file supplied by 'genome_file' parameter. Attempting automated download from MicrobesOnline...")

        # download genome from microbes online. store in MongoDB collection
        url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=genome" % ncbi_code
        save_name = "%d_genome.fa" % ncbi_code
        _download_url(url, save_name)

        seqs_b = []
        with open(save_name, 'r') as f:
            fasta_sequences = SeqIO.parse(f ,'fasta')
            for fasta in fasta_sequences:
                seqs_b.append({"scaffoldId": fasta.id,
                               "NCBI_RefSeq": fasta.description.split(" ")[1],
                               "NCBI_taxonomyId": fasta.description.split(" ")[-1],
                               "sequence": str(fasta.seq)})
    else:
        seqs_b = []

    # Check whether documents are already present in the collection before insertion
    if db.genome.count() > 0:
        seqs_f = filter(lambda doc: _not_exists(db.genome, doc), seqs_b)
    else:
        seqs_f = seqs_b

    logging.info("%d new genome sequence records to write", len(seqs_f))
    if len(seqs_f) > 0:
        db.genome.insert(seqs_f)

    return db.genome


def _load_ratios(file_in):
    """Loads ratios from individual cMonkey runs (unfinished) or single flat file (gzip compressed)."""
    if file_in == None:
        # compile from individual runs
        # do be done
        file_in = sorted(glob.glob(os.path.join(ensembledir, "%s???/ratios.tsv.gz" % prefix)))
    else:
        logging.info("Loading gene expression file from '%s'", file_in)
        # load directly from gzip file
        ratios = pd.read_csv(gzip.open(file_in, 'rb'), index_col=0, sep="\t")

        if ratios.shape[1] == 0:  # attempt using comma as a delimiter if tab failed
             ratios = pd.read_csv(gzip.open(file_in, 'rb'), index_col=0, sep=",")

        if ratios.shape[1] == 0:
            # still wrong delimiter
            raise Exception("Cannot read ratios file. Check delimiter. Should be '\t' or ',' ")
    return ratios

def _standardize_ratios(ratios):
    """compute standardized ratios (global). row standardized"""
    ratios_standardized = ratios.copy()
    zscore = lambda x: (x - x.mean()) / x.std()

    for row in ratios.iterrows():
        ratios_standardized.loc[row[0]] = zscore(row[1])

    return ratios_standardized


def _cond_info2dict(col_table, cond_name):
    cond_data = col_table[col_table.egrin2_col_name == cond_name]
    cond_dict = {"col_id": np.unique(cond_data.col_id)[0], "egrin2_col_name": np.unique(cond_data.egrin2_col_name)[0], "additional_info": []}

    if "feature_name" in cond_data.columns and "value" in cond_data.columns and "feature_units" in cond_data.columns:
        for i in range(0, cond_data.shape[0]):
            cond_dict["additional_info"].append({"name": cond_data.irow(i)["feature_name"],
                                                 "value": cond_data.irow(i)["value"],
                                                 "units": cond_data.irow(i)["feature_units"]})
    return cond_dict



###########################################################################
### Result Database class
### Coordinates the collection and assembly of corem information
###########################################################################

class ResultDatabase:

    def __init__(self, organism, db, ensembledir, prefix, ratios_raw, gre2motif, col_annot,
                 ncbi_code, genome_file, row_annot, row_annot_match_col,
                 targetdir, db_run_override=False):
        self.organism = organism
        self.db = db
        self.prefix = prefix
        self.ensembledir = ensembledir
        self.gre2motif = _gre2motif_path(ensembledir, gre2motif)
        self.db_files = _available_cmonkey_resultdb_files(ensembledir, prefix)
        self.db_run_override = db_run_override
        self.ncbi_code = _get_ncbi_code(ncbi_code, self.db_files[0])
        self.ratios_raw = ratios_raw
        self.col_annot = col_annot
        self.genome_file = genome_file
        self.row_annot = row_annot
        self.row_annot_match_col = row_annot_match_col

    def __get_row2id(self, ratios_standardized):
        """make row2id and id2row dicts for lookup"""
        row_name = []
        row_id = []

        for i in self.db.row_info.find():
            row_name.append(i["egrin2_row_name"])
            row_id.append(i["row_id"])

        for i in ratios_standardized.index.values:
            if i not in row_name:
                row_name.append(i)
            if len(row_id) > 0:
                row_id.append(max(row_id) + 1)
            else:
                row_id.append(0)

        row_info =  pd.DataFrame(zip(row_id, row_name), index=row_name, columns=["row_id", "egrin2_row_name"])
        return row_info

    def __insert_row_info(self):
        """
        Insert row_info into mongoDB database

        example queries
        ------------------------------
        for i in row_info_collection.find( { "name" : "carA" } ):
            print i
        for i in row_info_collection.find( { "sysName": { "$in" : ["b0032", "b0124", "b0089", "b0432","b2234","b0456"] } } ):
            print "%s = %s" % ( i["sysName"], i["name"] )
            print "It does: %s" % i["GO"]

        for i in row_info_collection.find( { "GO": { "$regex" : "GO:0006541," } } ):
            print i["sysName"]
            print i["GO"]

        """
        row_annot = self.row_annot
        row_annot_match_col = self.row_annot_match_col

        if row_annot is None:
            logging.info("Row annotation file not supplied as 'row_annot' parameter. Attempting automated download from MicrobesOnline...")
            if self.ncbi_code is not None:
                logging.info("Downloading gene information for NCBI taxonomy ID: %s", str(self.ncbi_code))
                url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%s;export=tab" % str(self.ncbi_code)
                save_name = "%s_geneInfo.tab" % str(self.ncbi_code)
                _download_url(url, save_name)

                with open(save_name, 'r') as f:
                    row_annot = pd.read_csv(f, sep="\t")

                left_on = "egrin2_row_name"

                if row_annot_match_col is None:
                    row_annot_match_col = "sysName"

                # join with row_annot
                row_table = pd.merge(self.row2id, row_annot, left_on=left_on, right_on=row_annot_match_col)
            else:
                print "WARNING: could not fetch additional gene information for NCBI taxonomy ID:", self.ncbi_code
                # TODO:
                # Check whether documents are already present in the collection before insertion
                # In case where they are present, update them
                row_table = self.row2id
        else:
            row_annot = pd.read_csv(open(row_annot, 'rb'), sep="\t")
            # join with row_annot
            row_table = pd.merge(self.row2id, row_annot, left_on=left_on, right_on=row_annot_match_col)

        # Check whether documents are already present in the collection before insertion
        docs = row_table.to_dict('records')

        if self.db.row_info.count() > 0:
            d_f = filter(lambda doc: _not_exists(self.db.row_info, doc), docs)
        else:
            d_f = docs

        logging.info("%d new row info records to write", len(d_f))

        if len(d_f) > 0:
            self.db.row_info.insert(d_f)

        return self.db.row_info

    def __get_col2id(self, ratios_standardized):
        """make cond2id and id2cond dicts for lookup"""
        col_name = []
        col_id = []

        for i in self.db.col_info.find():
            col_name.append(i["egrin2_col_name"])
            col_id.append(i["col_id"])

        for i in ratios_standardized.columns.values:
            if i not in col_name:
                col_name.append(i)
                if len(col_id) > 0:
                    col_id.append(max(col_id) + 1)
                else:
                    col_id.append(0)
        col_info = pd.DataFrame(zip(col_id, col_name), index=col_name, columns=["col_id", "egrin2_col_name"])
        return col_info

    def __insert_col_info(self):
        """
        Insert col_info into mongoDB database

        example queries
        ------------------------------
        for i in col_info_collection.find( { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ):
            print i["egrin2_col_name"]

        for i in col_info_collection.find( { "$and": [ { "additional_info.name": "strain", "additional_info.value":  { "$regex": "MG1655" } }, { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ] } ):
            print i
        """
        col_annot = self.col_annot

        # load additional col info
        if col_annot is not None:
            # assumes default microbes online schema
            col_annot = pd.read_csv(gzip.open(col_annot, 'rb'), sep="\t")
            col_table = pd.merge(self.col2id, col_annot, left_on="egrin2_col_name", right_on="experiment_name")
        else:
            col_table = self.col2id

        col_info_4_mongoDB = []
        for i in range(0, len(self.col2id.egrin2_col_name)):
            col_info_4_mongoDB.append(_cond_info2dict(col_table, self.col2id.egrin2_col_name[i]))

        # Check whether documents are already present in the collection before insertion
        if self.db.col_info.count() > 0:
            d_f = filter(lambda doc: _not_exists(self.db.col_info, doc), col_info_4_mongoDB)
        else:
            d_f = col_info_4_mongoDB

        logging.info("%d new col info records to write", len(d_f))

        if len(d_f) > 0:
            self.db.col_info.insert(d_f)

        return self.db.col_info

    def __insert_gene_expression(self):
        """
        Insert gene_expression into mongoDB database
        """
        ratios = self.expression
        ratios_std = self.expression_standardized

        exp_data = []
        num_rows = 0

        # speed up the lookups by putting this into a dictionary by replacing row/column
        # lookups with dictionaries and pulling out per-row lookups out of the inner loop

        # TODO: Pandas dataframes are probably not a good idea for lookup tables
        # ----- get rid of them for row2id/col2id altogether
        row_map = {row: self.row2id.loc[row].row_id for row in ratios.index.values}
        col_map = {col: self.col2id.loc[col].col_id for col in ratios.columns.values}
        
        for row_name in ratios.index.values:
            if num_rows % 200 == 0:
                logging.info("%.2f percent done (%d rows)", round((float(num_rows) / ratios.shape[0]) * 100, 1), num_rows)
            
            raw_row = ratios.loc[row_name]
            std_row = ratios_std.loc[row_name]
            row_id = row_map[row_name]

            for col_name in ratios.columns.values:
                raw_exp = raw_row[col_name]
                std_exp = std_row[col_name]
                exp_data.append({"row_id": row_id,
                                 "col_id": col_map[col_name],
                                 "raw_expression": raw_exp,
                                 "standardized_expression": std_exp
                                 })
            num_rows += 1

        # write to mongoDB collection
        gene_expression_collection = self.db.gene_expression

        # Check whether documents are already present in the collection before insertion
        if self.db.gene_expression.count() > 0:
            d_f = filter(lambda doc: _not_exists(self.db.gene_expression, doc), exp_data)
        else:
            d_f = exp_data

        logging.info("%d new gene expression records to write", len(d_f))

        if len(d_f) > 0:
            self.db.gene_expression.insert(d_f)

        return self.db.gene_expression

    def __assemble_ensemble_info(self, db_file):
        """Create python ensemble_info dictionary for bulk import into MongoDB collections"""
        run_name = db_file.split("/")[-2]
        logging.info("Assembling run info for cMonkey run: %s", run_name)
        conn = sqlite3.connect(db_file)
        try:
            c = conn.cursor()
            c.execute("SELECT start_time, finish_time, num_iterations, organism, species, num_rows, num_columns, num_clusters, git_sha FROM run_infos")
            run_info = c.fetchone()
            c.execute("SELECT name FROM row_names")
            rows = [self.row2id.loc[str(row[0])].row_id for row in c.fetchall()]
            c.execute("SELECT name FROM column_names")
            cols = [self.col2id.loc[str(row[0])].col_id for row in c.fetchall()]

            d = {
                "run_id": self.run2id.loc[run_name].run_id,
                "run_name": run_name,
                "start_time": str(run_info[0]),
                "finish_time": str(run_info[1]),
                "num_iterations": int(run_info[2]),
                "organism": str(run_info[3]),
                "species": str(run_info[4]),
                "num_rows": int(run_info[5]),
                "rows": rows,
                "num_columns": int(run_info[6]),
                "cols": cols,
                "num_clusters": int(run_info[7]),
                "git_sha": str(run_info[8]),
                "added_to_ensemble": datetime.datetime.utcnow()
            }
        finally:
            conn.close()
        return d

    def __insert_ensemble_info(self):
        """Compile and insert ensemble_info collection into MongoDB collection"""
        to_insert = [self.__assemble_ensemble_info(dbfile) for dbfile in self.db_files]
        ensemble_info_collection = self.db.ensemble_info

        # Check whether documents are already present in the collection before insertion
        if ensemble_info_collection.count() > 0:
            d_f = filter(lambda doc: _not_exists_key_value(ensemble_info_collection,
                                                           "run_name", doc["run_name"]), to_insert)
        else:
            d_f = to_insert

        logging.info("%d new ensemble info records to write", len(d_f))

        if len(d_f) > 0:
            ensemble_info_collection.insert(d_f)

        return ensemble_info_collection

    def __load_gre_map(self):
        if self.gre2motif is not None:
            count = 1
            mots = {}
            with open(self.gre2motif, 'r') as f:
                for line in f:
                    # only consider motif clusters with > 3 motifs
                    if len(line.strip("\n").split("\t")) > 3:
                        for motif in line.strip("\n").split("\t"):
                            elements = motif.split("_")
                            # cluster
                            elements[1] = int(elements[1])
                            # motif_num
                            elements[2] = int(elements[2])
                            # elements[0] = run_name
                            if elements[0] in mots.keys():
                                if elements[1] in mots[elements[0]].keys():
                                    mots[elements[0]][elements[1]][elements[2]] = count
                                else:
                                    mots[elements[0]][elements[1]] = {}
                                    mots[elements[0]][elements[1]][elements[2]] = count
                            else:
                                mots[elements[0]] = {}
                                mots[elements[0]][elements[1]] = {}
                                mots[elements[0]][elements[1]][elements[2]] = count
                        count += 1
            return mots
        else:
            return None

    def __insert_bicluster_info(self, db_file):
        """Find all biclusters in a cMonkey run, process and add as documents to bicluster collection

        example queries
        ------------------------------
        bicluster_info_collection.find({"rows":{"$all":[26,27]}}).count()
        """
        # Get all biclusters from cmonkey run
        conn = sqlite3.connect(db_file)
        try:
            c = conn.cursor()
            c.execute("SELECT max(iteration) FROM cluster_stats")
            last_run = c.fetchone()[0] # i think there is an indexing problem in cMonkey python!!
            w = (last_run, )
            c.execute("SELECT cluster FROM cluster_stats WHERE iteration = ?", w)
            biclusters = [self.__assemble_bicluster_info_single(db_file, c, last_run, row[0])
                          for row in c.fetchall()]
        finally:
            conn.close()

        bicluster_info_collection = self.db.bicluster_info

        # Check whether documents are already present in the collection before insertion
        if bicluster_info_collection.count() > 0:
            d_f = filter(lambda doc: _not_exists_key_value2(bicluster_info_collection,
                                                            "run_id", doc["run_id"],
                                                            "cluster", doc["cluster"]), biclusters)
        else:
            d_f = biclusters

        logging.info("%d new bicluster info records to write", len(d_f))
        if len(d_f) > 0:
            bicluster_info_collection.insert(d_f)

        return bicluster_info_collection

    def __assemble_bicluster_info_single(self, db_file, cursor, iteration, cluster):
        """Create python ensemble_info dictionary for bulk import into MongoDB collections"""
        run_name = db_file.split("/")[-2]
        w = (cluster, iteration)
        cursor.execute("SELECT residual FROM cluster_stats WHERE cluster = ? AND iteration = ?", w)
        residual = cursor.fetchone()[0]
        cursor.execute("SELECT name FROM row_members JOIN row_names ON row_members.order_num = row_names.order_num WHERE row_members.cluster = ? AND row_members.iteration = ?", w)
        rows = [self.row2id.loc[str(i[0])].row_id for i in cursor.fetchall()]
        cursor.execute("SELECT name FROM column_members JOIN column_names ON column_members.order_num = column_names.order_num WHERE column_members.cluster = ? AND column_members.iteration = ?", w)
        cols = [self.col2id.loc[str(i[0])].col_id for i in cursor.fetchall()]

        d = {
            "run_id": self.run2id.loc[run_name].run_id,
            "cluster": cluster,
            "rows": rows,
            "columns": cols,
            "residual": residual,
        }

        return d

    def __insert_motif_info(self, db_file):
        # Get all biclusters from cmonkey run
        conn = sqlite3.connect(db_file)
        try:
            c = conn.cursor()
            c.execute("SELECT max(iteration) FROM cluster_stats")
            last_run = c.fetchone()[0] # i think there is an indexing problem in cMonkey python!!
            w = (last_run, )

            c.execute("SELECT cluster FROM cluster_stats WHERE iteration=?", w)
            d_f = [self.__assemble_motif_info_single(db_file, c, last_run, row[0])
                   for row in c.fetchall()]
        finally:
            conn.close()

        d_f = list(itertools.chain(*d_f))
        d_f = filter(None, d_f)

        if len(d_f) > 0:
            self.db.motif_info.insert(d_f)

        return None

    def __assemble_motif_info_single(self, db_file, cursor, iteration, cluster):
        run_name = db_file.split("/")[-2]
        cluster_id = list(self.db.bicluster_info.find({"run_id": self.run2id.loc[run_name].run_id, "cluster": cluster}, {"_id": 1}))
        if len(cluster_id) > 1:
            logging.info("Cluster %s from run %s matches more than one entry in MongoDB. You have a problem. Refusing to add motif_info.", cluster, run_name)
            return None
        else:
            cluster_id = ObjectId(cluster_id[0]["_id"])

        cursor.execute("select motif_num from motif_infos where cluster=? and iteration=?", [cluster, iteration])
        motif_nums = [row[0] for row in cursor.fetchall()]
        motif_info = [self.__get_motif_info_single(cursor, iteration, run_name, cluster, i, cluster_id)
                      for i in motif_nums]

        return motif_info

    def __get_motif_info_single(self, cursor, iteration, run_name, cluster, motif_num, cluster_id):
        w = (cluster, iteration, motif_num)
        cursor.execute("SELECT seqtype, evalue FROM motif_infos WHERE cluster = ? AND iteration = ? AND motif_num = ?", w)
        motif_data = cursor.fetchone()
        cursor.execute("SELECT meme_motif_sites.rowid FROM meme_motif_sites JOIN motif_infos ON meme_motif_sites.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ?", w)
        meme_site = [i[0] for i in cursor.fetchall()]
        cursor.execute("SELECT row FROM motif_pssm_rows JOIN motif_infos ON motif_pssm_rows.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ?", w)
        rows = [i[0] for i in cursor.fetchall()]

        try:
            gre_id = self.motif2gre[run_name][cluster][motif_num]
        except:
            gre_id = "NaN"

        d = {
            "cluster_id": cluster_id,
            "gre_id": gre_id,
            "motif_num": motif_num,
            "seqtype": motif_data[0],
            "evalue": motif_data[1],
            "meme_motif_site": [self.__get_meme_motif_site_single(cursor, iteration, run_name, cluster, motif_num, i)
                                                                  for i in meme_site],
            "pwm": [self.__get_pwm_single(cursor, iteration, run_name, cluster, motif_num, i)
                    for i in rows]
        }
        return d

    def __get_meme_motif_site_single(self, cursor, iteration, run_name, cluster, motif_num, rowid):
        w = (str(rowid), )
        cursor.execute("SELECT seq_name, reverse, start, pvalue, flank_left, seq, flank_right FROM meme_motif_sites WHERE rowid = ?", w)
        data = cursor.fetchone()

        # try to match accession to row_id
        try:
            # translate from accession to row_id, requires microbes online
            row_id = self.row_info_collection.find({"accession": data[0]})[0]["row_id"]
        except:
            row_id = "NaN"

        try:
            scaffoldId = self.row_info_collection.find({"accession": data[0]})[0]["scaffoldId"]
        except:
            scaffoldId = "NaN"

        d = {
            "row_id": row_id,
            "reverse": data[1],
            "scaffoldId": scaffoldId,
            "start": data[2],
            "pvalue": data[3],
        }
        return d

    def __get_pwm_single(self, cursor, iteration, run_name, cluster, motif_num, row):
        w = (cluster, iteration, motif_num, row)
        cursor.execute("SELECT row, a, c, g, t FROM motif_pssm_rows JOIN motif_infos ON motif_pssm_rows.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ? AND motif_pssm_rows.row = ?", w)
        data = cursor.fetchone()

        # TODO: WW: This format for PSSMs is not really efficient for extraction
        # ---- We should rewrite it into arrays
        d = {
            "row": data[0],
            "a": data[1],
            "c": data[2],
            "g": data[3],
            "t": data[4]
        }
        return d

    def __assemble_fimo(self):

        def get_fimo_scans_single(i, db, ensembledir, run2id):
            cluster = i.cluster
            run_name =run2id[run2id["run_id"] == i.run_id]["run_name"][0]
            cluster_id = i._id

            try:
                # get all fimo scans in the dir
                f = os.path.join(ensembledir, run_name, "fimo-outs", "fimo-out-%04d.bz2" % cluster)
                fimo = pd.read_csv(f, sep="\t", compression = "bz2")

                # change sequence_name to scaffoldId
                fimo.rename(columns={'matched sequence': 'matched_sequence', 'sequence name': 'scaffoldId', "#pattern name": "motif_num"}, inplace=True)

                trans_d = {}
                for i in np.unique(fimo.scaffoldId):
                    NCBI_RefSeq = "_".join(i.split('.')[-2].split("_")[:: -1][0: 2][:: -1])
                    scaffoldId = db.genome.find_one({"NCBI_RefSeq": NCBI_RefSeq})["scaffoldId"]
                    trans_d[i] = scaffoldId

                trans_v = [trans_d[i] for i in fimo.scaffoldId.values]
                fimo.scaffoldId = trans_v
                fimo["cluster_id"] = cluster_id

                # only keep specific columns
                fimo = fimo.loc[:, ['scaffoldId', 'start', 'stop', 'strand', 'score', 'p-value', 'in_coding_rgn', 'cluster_id', 'motif_num']]

                # only store hits with p-value lte 1e-5
                fimo = fimo.loc[fimo["p-value"] <= 1e-5, ]
                d_f = fimo.to_dict('records')
                db.fimo.insert(d_f)

                # insert into fimo_small only if the motif maps to a GRE and the pval is less than 1e-5
                mot2gre = pd.DataFrame(list(db.motif_info.find({"cluster_id": cluster_id}, {"motif_num": 1, "gre_id": 1})))

                for x in range( mot2gre.shape[0]):
                    if mot2gre.iloc[x].gre_id != "NaN":
                        tmp_fimo = fimo.loc[fimo["motif_num"] == mot2gre.iloc[x].motif_num, ]
                        d_f = tmp_fimo.to_dict('records')
                        db.fimo_small.insert(d_f)
                return None
            except Exception:
                return None

        # get all biclusters
        bcs = pd.DataFrame(list( self.db.bicluster_info.find({}, {"cluster": 1, "run_id": 1})))
        tmp = bcs.apply(get_fimo_scans_single, axis=1, db=self.db, ensembledir=self.ensembledir, run2id=self.run2id)
        return None

    def mongo_dump(self, db, outfile, add_files=None):
        """Write contents from MongoDB instance to binary file"""
        logging.info("Dumping MongoDB to BSON")
        outfile_wdir = os.path.abspath(os.path.join(self.targetdir, outfile))
        self.db.client
        sys_command = "mongodump --db %s --out %s" % (db, outfile_wdir)
        os.system(sys_command)

        logging.info("Compressing MongoDB BSON docs")
        if add_files is not None:
            sys_command2 = "tar -czvf " + outfile + ".tgz " + "-C " + os.path.abspath(self.targetdir) + " " + outfile + " " + add_files
        else:
            sys_command2 = "tar -czvf " + outfile + ".tgz " + "-C " + os.path.abspath(self.targetdir) + " " + outfile
        os.system(sys_command2)

        logging.info("Cleaning up...")
        sys_command3 = "rm -rf " + outfile_wdir
        os.system(sys_command3)

    def mongo_restore(self, db, infile):
        """Read contents of binary MongoDB dump into MongoDB instance"""
        sys_command = "mongorestore --db %s %s" % (db, infile)
        os.system(sys_command)

    def compile(self):
        """Compile EGRIN2 ensemble"""
        logging.info("Compiling EGRIN2 ensemble...")
        self.db_files = _valid_cmonkey_results(self.db, self.db_files, self.db_run_override)
        self.run2id = _get_run2id(self.db, self.db_files)

        self.genome_collection = _import_genome(self.db, self.genome_file, self.ncbi_code)
        self.expression = _load_ratios(self.ratios_raw)

        logging.info("Standardizing gene expression...")
        self.expression_standardized = _standardize_ratios(self.expression)

        logging.info("Inserting into row_info collection")
        self.row2id = self.__get_row2id(self.expression_standardized)
        self.row_info_collection = self.__insert_row_info()

        logging.info("Inserting into col_info collection")
        self.col2id = self.__get_col2id(self.expression_standardized)
        self.col_info_collection = self.__insert_col_info()

        logging.info("Inserting gene expression into database")
        self.gene_expression_collection = self.__insert_gene_expression()
        self.gene_expression_collection.ensure_index("row_id")
        self.gene_expression_collection.ensure_index("col_id")

        logging.info("Inserting into ensemble_info collection")
        self.ensemble_info_collection = self.__insert_ensemble_info()
        self.motif2gre = self.__load_gre_map()

        logging.info("Inserting into bicluster collection")
        for dbfile in self.db_files:
            logging.info("%s", str(dbfile))
            self.bicluster_info_collection = self.__insert_bicluster_info(dbfile)
            self.__insert_motif_info(dbfile)

        logging.info("Indexing bicluster collection")
        self.bicluster_info_collection.ensure_index("rows")
        self.bicluster_info_collection.ensure_index("columns")
        self.db.motif_info.ensure_index("cluster_id")
        self.db.motif_info.ensure_index("gre_id")

        logging.info("Inserting into fimo collection. This might take awhile...")
        self.__assemble_fimo()

        logging.info("Indexing fimo collection")
        self.db.fimo.ensure_index([("scaffoldId", pymongo.ASCENDING), ("start", pymongo.ASCENDING), ("stop", pymongo.ASCENDING),
                                   ("p-value", pymongo.ASCENDING), ("cluster_id", pymongo.ASCENDING)])
        self.db.fimo_small.ensure_index([("scaffoldId", pymongo.ASCENDING), ("start", pymongo.ASCENDING),
                                         ("stop", pymongo.ASCENDING), ("p-value", pymongo.ASCENDING),
                                         ("cluster_id", pymongo.ASCENDING)])
        return None
