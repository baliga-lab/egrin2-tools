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
from pymongo import MongoClient
from bson.objectid import ObjectId
from Bio import SeqIO


# ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz"
# col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz"
# db_file = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/eco-out-001/cmonkey_run.db"

# delete egrin2 db. run on console

# sql2mongoDB( ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz",  col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz", ncbi_code = "511145")

class ResultDatabase:

    def __init__(self, organism, host, port, ensembledir=None, targetdir=None,
                 prefix=None, ratios_raw=None, gre2motif=None, col_annot=None, ncbi_code=None,
                 dbname=None , db_run_override=None, genome_file=None, row_annot=None,
                 row_annot_match_col=None):

        self.organism = organism
        self.host = host
        self.port = port

        client = MongoClient(host=host, port=port)
        logging.info("Connected to MongoDB")
        
        self.dbname = '%s_db' % self.organism if dbname is None else dbname

        if self.dbname in client.database_names():
            logging.warn("WARNING: %s database already exists!!!", self.dbname)
        else:
            logging.info("Initializing MongoDB database: %s", self.dbname)

        self.db = client[self.dbname]

        # get db files in directory
        if prefix is None:
            self.prefix = organism + '-out-'
        else:
            self.prefix = prefix
        if ensembledir is None:
            self.ensembledir = '.'
        else:
            self.ensembledir = ensembledir
        if targetdir is None:
            self.targetdir = '.'
        else:
            self.targetdir = targetdir

        if gre2motif == None:
            # default file name?
            if os.path.isfile(os.path.join(self.ensembledir, "out.mot_metaclustering.txt.I45.txt")):
                self.gre2motif = os.path.join(self.ensembledir, "out.mot_metaclustering.txt.I45.txt")
            else:
                logging.warn("I cannot find a GRE clustering file. If you want to assign GREs, please specify this file.")
                self.gre2motif = None
        else:
            self.gre2motif = gre2motif

        # get all cmonkey_run.db files
        self.db_files = np.sort(np.array(glob.glob(os.path.join(self.ensembledir, self.prefix) + "???/cmonkey_run.db")))
        self.db_run_override = db_run_override

        if ncbi_code == None:
            # try to find ncbi_code in cmonkey_run.db
            conn = sqlite3.connect(self.db_files[0])
            c = conn.cursor()
            try:
                c.execute("SELECT ncbi_code FROM run_infos")
                self.ncbi_code = c.fetchone()[0]
            except sqlite3.Error as e:
                logging.warn("Could not find NCBI Genome ID in cmonkey_run.db: %s", e.args[0])
        else:
            self.ncbi_code = ncbi_code

        self.ratios_raw = ratios_raw
        self.col_annot = col_annot
        self.genome_file = genome_file
        self.row_annot = row_annot
        self.row_annot_match_col = row_annot_match_col

        if len(self.db_files) < 1:
            logging.warn("""I cannot find any cMonkey SQLite databases in the current directory: %s
Make sure 'ensembledir' variable points to the location of your cMonkey-2 ensemble results.""",
                         os.getcwd())
        return None

    def dlfile(self, url, save_name=None, num_retries=5):
        # Open the url
        if save_name is None:
            save_name = "tmp"

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

    def checkRuns(self, db_files, db_run_override, db):
        """make sure the runs have data!!!"""
        to_keep = []
        ensemble_info_collection = db.ensemble_info
        for i in db_files:
            try:
                conn = sqlite3.connect(i)
                c = conn.cursor()
                c.execute("SELECT num_iterations FROM run_infos")
                run_info = c.fetchone()

                if run_info is None:
                    pass
                elif run_info[0] < 1000:
                    # incomplete runs
                    pass
                else:
                    if db_run_override == None:
                        # do not include runs that are already in the database
                        # check for existence
                        run_name = i.split("/")[-2]
                        if ensemble_info_collection.find({"run_name": run_name}).count() > 0:
                            pass
                        else:
                            to_keep.append(i)
                    else:
                        to_keep.append(i)
            except Exception:
                pass
        return to_keep

    def get_run2id(self, dbfiles, db):
        """make run2id"""
        ensemble_info_collection = db.ensemble_info
        run_name = []
        run_id = []

        for i in ensemble_info_collection.find():
            run_name.append(i["run_name"])
            run_id.append(i["run_id"])

        for i in dbfiles:
            if i.split("/")[-2] not in run_name:
                run_name.append(i.split("/")[-2])
                if len(run_id) > 0:
                    run_id.append(max(run_id) + 1)
                else:
                    run_id.append(0)
        run_info =  pd.DataFrame(zip(run_id, run_name), index=run_name, columns=["run_id", "run_name"])
        return run_info

    def check4existence(self, collection, document, key1=None, value1=None, key2=None, value2=None):
        if key1 == None:
            d_check = collection.find(document).count()
        else:
            if key2 == None:
                d_check = collection.find({key1: value1}).count()
            else:
                d_check = collection.find({key1: value1, key2: value2}).count()
        if d_check == 0:
            return document

    def loadGenome (self, ncbi_code, genome_file):
        if genome_file == None:
            logging.info("No custom genome annotation file supplied by 'genome_file' parameter. Attempting automated download from MicrobesOnline...")

            # download genome from microbes online. store in MongoDB collection
            url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%s;export=genome" % str(ncbi_code)
            save_name = "%s_genome.fa" % str(ncbi_code)
            self.dlfile(url, save_name)

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

        genome_collection = self.db.genome

        # Check whether documents are already present in the collection before insertion
        if genome_collection.count() > 0:
            seqs_f = filter(None, [self.check4existence(genome_collection, i) for i in seqs_b])
        else:
            seqs_f = seqs_b

        logging.info("%d new genome sequence records to write", len(seqs_f))
        if len(seqs_f) > 0:
            genome_collection.insert(seqs_f)

        return genome_collection

    def loadRatios(self, file_in):
        """Loads ratios from individual cMonkey runs (unfinished) or single flat file (gzip compressed)."""
        if file_in == None:
            # compile from individual runs
            # do be done
            file_in = np.sort(np.array(glob.glob(ensembledir + prefix + "???/ratios.tsv.gz")))
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

    def standardizeRatios(self, ratios):
        """compute standardized ratios (global). row standardized"""
        ratios_standardized = ratios.copy()
        zscore = lambda x: (x - x.mean()) / x.std()

        for row in ratios.iterrows():
            ratios_standardized.loc[row[0]] = zscore(row[1])

        return ratios_standardized

    def get_row2id(self, ratios_standardized, db):
        """make row2id and id2row dicts for lookup"""
        row_info_collection = db.row_info
        row_name = []
        row_id = []

        for i in row_info_collection.find():
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

    def insert_row_info(self, ncbi_code, row_info, row_annot, row_annot_match_col):
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
        if row_annot == None:
            logging.info("Row annotation file not supplied as 'row_annot' parameter. Attempting automated download from MicrobesOnline...")
            if ncbi_code is not None:
                logging.info("Downloading gene information for NCBI taxonomy ID: %s", str(ncbi_code))
                url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%s;export=tab" % str(ncbi_code)
                save_name = "%s_geneInfo.tab" % str(ncbi_code)
                self.dlfile(url, save_name)

                with open(save_name, 'r') as f:
                    row_annot = pd.read_csv(f, sep="\t")

                left_on = "egrin2_row_name"

                if row_annot_match_col == None:
                    row_annot_match_col = "sysName"

                # join with row_annot
                row_table = pd.merge(row_info, row_annot, left_on=left_on, right_on=row_annot_match_col)
            else:
                print "WARNING: could not fetch additional gene information for NCBI taxonomy ID:", ncbi_code
                # TODO:
                # Check whether documents are already present in the collection before insertion
                # In case where they are present, update them
                row_table = row_info
        else:
            row_annot = pd.read_csv(open(row_annot, 'rb'), sep="\t")
            # join with row_annot
            row_table = pd.merge(row_info, row_annot, left_on=left_on, right_on=row_annot_match_col)

        # write to mongoDB collection
        row_info_collection = self.db.row_info

        # Check whether documents are already present in the collection before insertion
        d = row_table.to_dict('records')

        if row_info_collection.count() > 0:
            d_f = filter(None, [self.check4existence( row_info_collection, i) for i in d])
        else:
            d_f = d

        logging.info("%d new row info records to write", len(d_f))

        if len(d_f) > 0:
            row_info_collection.insert(d_f)

        return row_info_collection

    def get_col2id(self, ratios_standardized, db):
        """make cond2id and id2cond dicts for lookup"""
        col_info_collection = db.col_info
        col_name = []
        col_id = []

        for i in col_info_collection.find():
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

    def insert_col_info(self, col_info, col_annot):
        """
        Insert col_info into mongoDB database


        example queries
        ------------------------------
        for i in col_info_collection.find( { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ):
            print i["egrin2_col_name"]

        for i in col_info_collection.find( { "$and": [ { "additional_info.name": "strain", "additional_info.value":  { "$regex": "MG1655" } }, { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ] } ):
            print i

        """
        # load additional col info
        if  col_annot != None:
            # assumes default microbes online schema
            col_annot = pd.read_csv( gzip.open(col_annot, 'rb'), sep="\t")
            col_table = pd.merge(col_info, col_annot, left_on="egrin2_col_name", right_on="experiment_name")
        else:
            col_table = col_info

        col_info_4_mongoDB = []
        for i in range(0, len(col_info.egrin2_col_name)):
            col_info_4_mongoDB.append(self.condInfo2Dict(col_table, col_info.egrin2_col_name[i]))

        # write to mongoDB collection
        col_info_collection = self.db.col_info

        # Check whether documents are already present in the collection before insertion
        if col_info_collection.count() > 0:
            d_f = filter(None, [self.check4existence( col_info_collection, i) for i in col_info_4_mongoDB])
        else:
            d_f = col_info_4_mongoDB

        logging.info("%d new col info records to write", len(d_f))

        if len(d_f) > 0:
            col_info_collection.insert(d_f)

        return col_info_collection

    def condInfo2Dict(self, col_table, cond_name):
        cond_data = col_table[col_table.egrin2_col_name == cond_name]
        cond_dict = {"col_id": np.unique(cond_data.col_id)[0], "egrin2_col_name": np.unique(cond_data.egrin2_col_name)[0], "additional_info": []}

        if "feature_name" in cond_data.columns and "value" in cond_data.columns and "feature_units" in cond_data.columns:
            for i in range(0, cond_data.shape[0]):
                cond_dict[ "additional_info" ].append({"name": cond_data.irow(i)["feature_name"],
                                                       "value": cond_data.irow(i)["value"],
                                                       "units": cond_data.irow(i)["feature_units"]})
        return cond_dict

    def insert_gene_expression(self, db, row2id, col2id, ratios, ratios_std):
        """
        Insert gene_expression into mongoDB database
        """
        exp_data = []
        num_rows = 0

        # speed up the lookups by putting this into a dictionary by replacing row/column
        # lookups with dictionaries and pulling out per-row lookups out of the inner loop

        # TODO: Pandas dataframes are probably not a good idea for lookup tables
        # ----- get rid of them for row2id/col2id altogether
        row_map = {row: row2id.loc[row].row_id for row in ratios.index.values}
        col_map = {col: col2id.loc[col].col_id for col in ratios.columns.values}
        
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
        gene_expression_collection = db.gene_expression

        # Check whether documents are already present in the collection before insertion
        if gene_expression_collection.count() > 0:
            d_f = filter(None, [self.check4existence(gene_expression_collection, i) for i in exp_data])
        else:
            d_f = exp_data

        logging.info("%d new gene expression records to write", len(d_f))

        if len(d_f) > 0:
            gene_expression_collection.insert(d_f)

        return gene_expression_collection

    def assemble_ensemble_info(self, db_file, run2id, row2id, col2id):
        """Create python ensemble_info dictionary for bulk import into MongoDB collections"""
        run_name = db_file.split("/")[-2]
        logging.info("Assembling run info for cMonkey run: %s", run_name)
        conn = sqlite3.connect(db_file)
        try:
            c = conn.cursor()
            c.execute("SELECT start_time, finish_time, num_iterations, organism, species, num_rows, num_columns, num_clusters, git_sha FROM run_infos")
            run_info = c.fetchone()
            c.execute("SELECT name FROM row_names")
            rows = [row2id.loc[str(row[0])].row_id for row in c.fetchall()]
            c.execute("SELECT name FROM column_names")
            cols = [col2id.loc[str(row[0])].col_id for row in c.fetchall()]

            d = {
                "run_id": run2id.loc[run_name].run_id,
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

    def insert_ensemble_info(self, db_files, db, run2id, row2id, col2id):
        """Compile and insert ensemble_info collection into MongoDB collection"""
        to_insert = [self.assemble_ensemble_info( i, run2id, row2id, col2id) for i in db_files]
        ensemble_info_collection = db.ensemble_info

        # Check whether documents are already present in the collection before insertion
        if ensemble_info_collection.count() > 0:
            d_f = filter(None, [self.check4existence(ensemble_info_collection, i, "run_name", i["run_name"])
                                for i in to_insert])
        else:
            d_f = to_insert

        logging.info("%d new ensemble info records to write", len(d_f))

        if len(d_f) > 0:
            ensemble_info_collection.insert(d_f)

        return ensemble_info_collection

    def loadGREMap(self, gre2motif):
        if gre2motif is not None:
            count = 1
            mots = {}
            with open(gre2motif, 'r') as f:
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

    def insert_bicluster_info( self, db, db_file, run2id, row2id, col2id ):
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
            biclusters = [self.assemble_bicluster_info_single(db, db_file, c, last_run, i[0], run2id, row2id, col2id)
                          for i in c.fetchall()]
        finally:
            conn.close()

        bicluster_info_collection = self.db.bicluster_info

        # Check whether documents are already present in the collection before insertion
        if bicluster_info_collection.count() > 0:
            d_f = filter(None, [self.check4existence(bicluster_info_collection, i, "run_id",
                                                     i["run_id"], "cluster", i["cluster"]) for i in biclusters])
        else:
            d_f = biclusters

        logging.info("%d new bicluster info records to write", len(d_f))
        if len(d_f) > 0:
            bicluster_info_collection.insert(d_f)

        return bicluster_info_collection

    def assemble_bicluster_info_single(self, db, db_file, cursor, iteration, cluster, run2id, row2id, col2id):
        """Create python ensemble_info dictionary for bulk import into MongoDB collections"""
        run_name = db_file.split("/")[-2]
        w = (cluster, iteration)
        cursor.execute("SELECT residual FROM cluster_stats WHERE cluster = ? AND iteration = ?", w)
        residual = cursor.fetchone()[0]
        cursor.execute("SELECT name FROM row_members JOIN row_names ON row_members.order_num = row_names.order_num WHERE row_members.cluster = ? AND row_members.iteration = ?", w)
        rows = [row2id.loc[str(i[0])].row_id for i in cursor.fetchall()]
        cursor.execute("SELECT name FROM column_members JOIN column_names ON column_members.order_num = column_names.order_num WHERE column_members.cluster = ? AND column_members.iteration = ?", w)
        cols = [col2id.loc[str(i[0])].col_id for i in cursor.fetchall()]

        d = {
            "run_id": run2id.loc[run_name].run_id,
            "cluster": cluster,
            "rows": rows,
            "columns": cols,
            "residual": residual,
        }

        return d

    def insert_motif_info(self, db, db_file, run2id, motif2gre, row_info_collection):
        # Get all biclusters from cmonkey run
        conn = sqlite3.connect(db_file)
        try:
            c = conn.cursor()
            c.execute("SELECT max(iteration) FROM cluster_stats")
            last_run = c.fetchone()[0] # i think there is an indexing problem in cMonkey python!!
            w = (last_run, )

            c.execute("SELECT cluster FROM cluster_stats WHERE iteration = ?", w)
            d_f = [self.assemble_motif_info_single(db, db_file, c, last_run, i[0], run2id, motif2gre, row_info_collection)
                   for i in c.fetchall()]
        finally:
            conn.close()

        d_f = list(itertools.chain(*d_f))
        d_f = filter(None, d_f)

        if len(d_f) > 0:
            db.motif_info.insert(d_f)

        return None

    def assemble_motif_info_single(self, db, db_file, cursor, iteration, cluster, run2id, motif2gre, row_info_collection):
        run_name = db_file.split("/")[-2]
        cluster_id = list(db.bicluster_info.find({"run_id": run2id.loc[run_name].run_id, "cluster": cluster}, {"_id": 1}))
        if len(cluster_id) > 1:
            logging.info("Cluster %s from run %s matches more than one entry in MongoDB. You have a problem. Refusing to add motif_info.", cluster, run_name)
            return None
        else:
            cluster_id = ObjectId(cluster_id[0]["_id"])

        w = (cluster, iteration)
        cursor.execute("SELECT motif_num FROM motif_infos WHERE cluster = ? AND iteration = ?", w)
        motif_nums = [i[0] for i in cursor.fetchall()]
        motif_info = [self.get_motif_info_single(db, cursor, iteration, run_name, cluster, i, motif2gre, row_info_collection, cluster_id) for i in motif_nums]

        return motif_info

    def get_motif_info_single(self, db, cursor, iteration, run_name, cluster, motif_num, motif2gre, row_info_collection, cluster_id):
        w = (cluster, iteration, motif_num)
        cursor.execute("SELECT seqtype, evalue FROM motif_infos WHERE cluster = ? AND iteration = ? AND motif_num = ?", w)
        motif_data = cursor.fetchone()
        cursor.execute("SELECT meme_motif_sites.rowid FROM meme_motif_sites JOIN motif_infos ON meme_motif_sites.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ?", w)
        meme_site = [i[0] for i in cursor.fetchall()]
        cursor.execute("SELECT row FROM motif_pssm_rows JOIN motif_infos ON motif_pssm_rows.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ?", w)
        rows = [i[0] for i in cursor.fetchall()]

        try:
            gre_id = motif2gre[run_name][cluster][motif_num]
        except:
            gre_id = "NaN"

        d = {
            "cluster_id": cluster_id,
            "gre_id": gre_id,
            "motif_num": motif_num,
            "seqtype": motif_data[0],
            "evalue": motif_data[1],
            "meme_motif_site": [self.get_meme_motif_site_single( cursor, iteration, run_name, cluster, motif_num, i, row_info_collection) for i in meme_site],
            "pwm": [self.get_pwm_single( cursor, iteration, run_name, cluster, motif_num, i, row_info_collection ) for i in rows]
        }
        return d

    def get_meme_motif_site_single(self, cursor, iteration, run_name, cluster, motif_num, rowid, row_info_collection):
        w = (str(rowid), )
        cursor.execute("SELECT seq_name, reverse, start, pvalue, flank_left, seq, flank_right FROM meme_motif_sites WHERE rowid = ?", w)
        data = cursor.fetchone()

        # try to match accession to row_id
        try:
            # translate from accession to row_id, requires microbes online
            row_id = row_info_collection.find({"accession": data[0]})[0]["row_id"]
        except:
            row_id = "NaN"

        try:
            scaffoldId = row_info_collection.find({"accession": data[0]})[0]["scaffoldId"]
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

    def get_pwm_single(self, cursor, iteration, run_name, cluster, motif_num, row, row_info_collection):
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

    def assemble_fimo(self):

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

    def mongoDump(self, db, outfile, add_files=None):
        """Write contents from MongoDB instance to binary file"""
        logging.info("Dumping MongoDB to BSON")
        outfile_wdir = os.path.abspath(os.path.join(self.targetdir, outfile))
        sys_command = "mongodump --db " + db + " --out " + outfile_wdir + " --host " + self.host + " --port " + str(self.port)
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

    def mongoRestore(self, db, infile):
        """Read contents of binary MongoDB dump into MongoDB instance"""
        sys_command = "mongorestore --db " + db + " --host " + self.host + " --port " + str(self.port) + " " + infile
        os.system(sys_command)

    def compile(self):
        """Compile EGRIN2 ensemble"""
        logging.info("Compiling EGRIN2 ensemble...")
        self.db_files = self.checkRuns(self.db_files, self.db_run_override, self.db)
        self.run2id = self.get_run2id(self.db_files, self.db)

        logging.info("Downloading genome information for NCBI taxonomy ID: %s", self.ncbi_code)
        self.genome_collection = self.loadGenome(self.ncbi_code, self.genome_file)
        self.expression = self.loadRatios(self.ratios_raw)

        logging.info("Standardizing gene expression...")
        self.expression_standardized = self.standardizeRatios(self.expression)

        logging.info("Inserting into row_info collection")
        self.row2id = self.get_row2id(self.expression_standardized, self.db)
        self.row_info_collection = self.insert_row_info(self.ncbi_code, self.row2id, self.row_annot, self.row_annot_match_col)

        logging.info("Inserting into col_info collection")
        self.col2id = self.get_col2id(self.expression_standardized, self.db)
        self.col_info_collection = self.insert_col_info(self.col2id, self.col_annot)

        logging.info("Inserting gene expression into database")
        self.gene_expression_collection = self.insert_gene_expression(self.db, self.row2id, self.col2id,
                                                                      self.expression, self.expression_standardized)
        self.gene_expression_collection.ensure_index("row_id")
        self.gene_expression_collection.ensure_index("col_id")

        logging.info("Inserting into ensemble_info collection")
        self.ensemble_info_collection = self.insert_ensemble_info(self.db_files, self.db, self.run2id, self.row2id, self.col2id)
        self.motif2gre = self.loadGREMap(self.gre2motif)

        logging.info("Inserting into bicluster collection")
        for i in self.db_files:
            logging.info("%s", str(i))
            self.bicluster_info_collection = self.insert_bicluster_info(self.db, i, self.run2id, self.row2id, self.col2id)
            self.insert_motif_info(self.db, i, self.run2id, self.motif2gre, self.row_info_collection)

        logging.info("Indexing bicluster collection")
        self.bicluster_info_collection.ensure_index("rows")
        self.bicluster_info_collection.ensure_index("columns")
        self.db.motif_info.ensure_index("cluster_id")
        self.db.motif_info.ensure_index("gre_id")

        logging.info("Inserting into fimo collection. This might take awhile...")
        self.assemble_fimo()

        logging.info("Indexing fimo collection")
        self.db.fimo.ensure_index([("scaffoldId", pymongo.ASCENDING), ("start", pymongo.ASCENDING), ("stop", pymongo.ASCENDING),
                                   ("p-value", pymongo.ASCENDING), ("cluster_id", pymongo.ASCENDING)])
        self.db.fimo_small.ensure_index([("scaffoldId", pymongo.ASCENDING), ("start", pymongo.ASCENDING),
                                         ("stop", pymongo.ASCENDING), ("p-value", pymongo.ASCENDING),
                                         ("cluster_id", pymongo.ASCENDING)])
        return None
