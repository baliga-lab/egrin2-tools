import os
import datetime
import glob
import sys
import gzip
import time
from urllib2 import urlopen, URLError, HTTPError
from zipfile import ZipFile

import pdb

import numpy as np
import pandas as pd

import sqlite3
import pymongo
from pymongo import MongoClient
import gridfs
from Bio import SeqIO

from egrin2.merge_egrin2_dbs import *

tmp = sql2mongoDB( ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz", col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz", ncbi_code = "511145",db_run_override = False)

# connect to database
# make sure mongodb is running
retvalue = os.system("nohup mongod --port 27017 --quiet &")
try:
	client = MongoClient('mongodb://localhost:27017/') 
	print "Connected to MongoDB"
except pymongo.errors.ConnectionFailure, e:
	print "Could not connect to MongoDB: %s" % e

dbname = "egrin2_db"

db = client[dbname]

# get db files in directory
prefix = 'eco-out-'
e_dir = './'

gre2motif = e_dir + "out.mot_metaclustering.txt.I45.txt"


db_files = tmp.db_files
db_files = tmp.checkRuns( db_files, tmp.db_run_override, db )
	
ncbi_code = "511145"

ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz"
col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz"

ratios = tmp.loadRatios( ratios_raw)
ratios_standardized =  tmp.standardizeRatios( ratios)

db_file = db_files[0]
run2id = tmp.get_run2id( db_files, db )
row2id = tmp.get_row2id( ratios_standardized, db )
col2id = tmp.get_col2id( ratios_standardized, db )
motif2gre = tmp.loadGREMap( gre2motif )
row_info_collection = db.row_info

conn = sqlite3.connect(db_file)
cursor = conn.cursor()

cursor.execute("SELECT max(iteration) FROM cluster_stats;")
iteration = cursor.fetchone()[0] # i think there is an indexing problem in cMonkey python!! 

cluster = 5

run_name = 'eco-out-001'

tmp.bicluster_info_collection = tmp.insert_bicluster_info( db, e_dir, db_file, run2id, row2id, col2id, motif2gre, row_info_collection )