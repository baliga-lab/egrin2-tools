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
import gridfs

# ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz"
# >>> import datetime
# >>> st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
# row_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/genomeInfo.tsv.gz"
# col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz"

class sql2mongoDB:
    
    def __init__( self, e_dir = None, prefix = None, gene_info = None, ratios_raw = None, col_annot = None, row_annot = None, dbname = None ):
    	
    	# connect to database
	# make sure mongodb is running
	# retvalue = os.system("nohup mongod --port 27017 --quiet &")
	try:
    		client = MongoClient('mongodb://localhost:27017/') 
    		print "Connected to MongoDB"
	except pymongo.errors.ConnectionFailure, e:
   		print "Could not connect to MongoDB: %s" % e

	if dbname = None:
		dbname = "egrin2_db"
	if dbname in client.database_names():
		print "WARNING: %s database already exists!!! You are probably just duplicating the data. I don't protect you from this currently." % dbname
	else:
		print "Initializing MongoDB database: %s" % dbname
	db = client[dbname]

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
        	self.ratios = loadRatios(self, ratios_raw)
        	self.ratios_standardized = standardizeRatios(self, self.ratios)
        	
        	#####################
        	# load additional row info
        	#####################
        	if  row_annot != None:
        		# assumes default microbes online schema
        		row_annot = pd.read_csv( gzip.open( row_annot, 'rb' ), sep="\t" )	
        	# make row2id and id2row dicts for lookup
        	rows = self.ratios_standardized.index.values
        	rows.sort()
        	self.row_info =  pd.DataFrame( zip( range( len( rows ) ), rows ), index = rows, columns = [ "row_id", "egrin2_row_name"] )
        	# join with row_annot
        	row_table = pd.merge( self.row_info, row_annot, left_on="egrin2_row_name", right_on="sysName" )
        	# write to mongoDB collection 
        	print "Inserting to row_info collection"
        	row_info_collection = db.row_info
        	row_info_collection.insert( row_table.to_dict( outtype='records' ) )
        	#
        	# example queries
        	#
        	# for i in row_info_collection.find( { "name" : "carA" } ):
        	# 	print i
        	# for i in row_info_collection.find( { "sysName": { "$in" : ["b0032", "b0124", "b0089", "b0432","b2234","b0456"] } } ):
        	# 	print "%s = %s" % ( i["sysName"], i["name"] )
        	# 	print "It does: %s" % i["GO"]
        	# for i in row_info_collection.find( { "GO": { "$regex" : "GO:0006541," } } ):
        	# 	print i["sysName"]
        	# 	print i["GO"]

        	#####################
        	# load additional col info
        	#####################
        	if  col_annot != None:
        		# assumes default microbes online schema
        		col_annot = pd.read_csv( gzip.open( col_annot, 'rb' ), sep="\t" )	
        	# make cond2id and id2cond dicts for lookup
        	print "Inserting to row_info collection"
        	cols = self.ratios_standardized.columns.values
        	cols.sort()
        	self.col_info =  pd.DataFrame( zip( range( len( cols ) ), cols ), columns = [ "col_id", "egrin2_col_name"] )
        	col_table = pd.merge( self.col_info, col_annot, left_on="egrin2_col_name", right_on="experiment_name" )
        	col_info_4_mongoDB = []
        	for i in range( 0,len( col_info.egrin2_col_name ) ):
        		col_info_4_mongoDB.append( condInfo2Dict( 0, col_table, col_info.egrin2_col_name[i] ) )
        	col_info_collection = db.col_info
        	col_info_collection.insert( col_info_4_mongoDB )
        	# for i in col_info_collection.find( { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ):
        	# 	print i["egrin2_col_name"]
        	# for i in col_info_collection.find( { "$and": [ { "additional_info.name": "strain", "additional_info.value":  { "$regex": "MG1655" } }, { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ] } ):
        	# 	print i
        	# for i in col_info_collection.find( { "$and": [ { "additional_info.value":   "0.0000303"  }, { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ] } ):
        	# 	print i["egrin2_col_name"]

    def loadRatios(self, file_in):
    	"""Loads ratios from individual cMonkey runs (unfinished) or single flat file (gzip compressed)."""
        	if file_in = None:
        		# compile from individual runs
        		# do be done
        		file_in = np.sort( np.array( glob.glob( e_dir + prefix + "???/ratios.tsv.gz" ) ) )
        	else:
        		print "Loading gene expression file from %s" % ratios_raw
        		# load directly from gzip file
        		ratios= pd.read_csv( gzip.open( file_in, 'rb' ), index_col=0, sep="\t" ) 
        	return ratios

     def standardizeRatios(self, ratios):
     	"""compute standardized ratios (global). row standardized"""
        	ratios_standardized = ratios.copy()
	zscore = lambda x: ( x - x.mean() ) / x.std()
	print "Standardizing gene expression..."
	for row in ratios.iterrows():
		ratios_standardized.loc[ row[0] ] = zscore( row[1] )
	return ratios_standardized
        	
    def condInfo2Dict( self, col_table, cond_name):
    	cond_data = col_table[col_table.egrin2_col_name==cond_name]
    	cond_dict = { "col_id": np.unique( cond_data.col_id )[0], "egrin2_col_name": np.unique( cond_data.egrin2_col_name )[0], "additional_info": [] }
    	for i in range(0,cond_data.shape[0]):
    		cond_dict[ "additional_info" ].append( { "name": cond_data.irow(i)["feature_name" ], "value": cond_data.irow(i)["value"], "units": cond_data.irow(i)["feature_units"] } )
    	return cond_dict

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




