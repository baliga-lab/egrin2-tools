#!/usr/bin/env python

"""Initialize  mongo database from individual cMonkey 
runs (SQLite) plus some additional tables and run_id column, 
including motifs and motif clusters"""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import os
import datetime
import glob
import sys
import gzip
import time

import pdb

import numpy as np
import pandas as pd

import sqlite3
from pymongo import MongoClient

# ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz"
# >>> import datetime
# >>> st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
# row_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/genomeInfo.tsv.gz"

class sql2mongoDB:
    
    def __init__( self, e_dir = None, prefix = None, gene_info = None, ratios_raw = None, col_annot = None, row_annot = None ):
    	
    	# connect to database
	# make sure mondod is running
	retvalue = os.system("nohup mongod --port 27017 --quiet &")

	client = MongoClient('mongodb://localhost:27017/') 

	db = client.egrin2_db

	collection = db.test_collection

    	# get db files in directory
    	if prefix is None:
    		prefix = 'eco-out-'
    	else:
        		prefix = prefix
        	if e_dir is None:
    		e_dir = './eco-ens-m3d/'
    	else:
        		prefix = prefix
        	self.db_files = np.sort( np.array( glob.glob( e_dir + prefix + "???/cmonkey_run.db" ) ) ) # get all cmonkey_run.db files

        	# make ratios
        	if ratios_raw = None:
        		# compile from individual runs
        		# do be done
        		ratios_files = np.sort( np.array( glob.glob( e_dir + prefix + "???/ratios.tsv.gz" ) ) )
        	else:
        		# load directly from gzip file
        		self.ratios_raw= pd.read_csv( gzip.open( ratios_raw, 'rb' ), index_col=0, sep="\t" ) 
        		self.ratios_standardized = self.ratios_raw.copy()
        		# compute standardized ratios (global)
        		# row standardized
        		zscore = lambda x: ( x - x.mean() ) / x.std()
        		for row in self.ratios_raw.iterrows():
        			self.ratios_standardized.loc[ row[0] ] = zscore( row[1] ) 

        	# load additional gene info
        	if  row_annot != None:
        		# assumes default microbes online schema
        		row_annot = pd.read_csv( gzip.open( row_annot, 'rb' ), index_col=7, sep="\t" )	
        	
        	# make row2id and id2row dicts for lookup
        	rows = self.ratios_standardized.index.values
        	rows.sort()
        	self.row_info =  pd.DataFrame( zip( range( len( rows ) ), rows ), index = rows, columns = [ "id", "egrin2_row_name"] )
        	#tmp.to_json(orient='records')

        	# make cond2id and id2cond dicts for lookup
        	cols = self.ratios_standardized.columns.values
        	cols.sort()
        	self.col_info =  pd.DataFrame( zip( range( len( cols ) ), cols ), columns = [ "id", "egrin2_col_name"] )
        	

    def parseRatios( self, ratios_files ):
    	"""Compile ratios from available runs. Note that each run stores only the ratios that were used for that run."""
    	## unnecessary at the moment since each of the standardized ratios are different at the moment
    	# i = 0
    	# while i < len(ratios_files):
    	# 	if i == 0:
    	# 		rats_df = pd.read_csv( gzip.open( ratios_files[i], 'rb' ), index_col=0, sep="\t" ) 
    	# 	else:
    	# 		new_df = pd.read_csv( gzip.open( ratios_files[i], 'rb' ), index_col=0, sep="\t" )
    	# 		rats_df = rats_df + new_df
    	# 	i = i+1 


    def parseSQLite( self, db_file ):
    	"""Create python dictionaries for bulk import into MongoDB collections"""
        	conn = sqlite3.connect(db_file)
        	c = conn.cursor()




