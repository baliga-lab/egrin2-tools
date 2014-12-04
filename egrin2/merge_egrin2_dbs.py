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
import pymongo
from pymongo im%pasteport MongoClient
import gridfs

# ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz"
# >>> import datetime
# >>> st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
# row_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/genomeInfo.tsv.gz"
# col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz"
# db_file = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/eco-out-001/cmonkey_run.db"

class sql2mongoDB:
    
	def __init__( self, e_dir = None, prefix = None, gene_info = None, ratios_raw = None, gre2motif = None, col_annot = None, row_annot = None, dbname = None ):
		
		# connect to database
		# make sure mongodb is running
		retvalue = os.system("nohup mongod --port 27017 --quiet &")
		try:
			client = MongoClient('mongodb://localhost:27017/') 
			print "Connected to MongoDB"
		except pymongo.errors.ConnectionFailure, e:
			print "Could not connect to MongoDB: %s" % e

		if dbname == None:
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
    		if gre2motif == None:
    			# default file name?
    			gre2motif = e_dir + "out.mot_metaclustering.txt.I45.txt"

	    	self.db_files = np.sort( np.array( glob.glob( e_dir + prefix + "???/cmonkey_run.db" ) ) ) # get all cmonkey_run.db files
	    	self.db_files = checkRuns( self, self.db_files )
	    	self.run2id = get_run2id(self, self.db_files )
	    	self.ratios = loadRatios(self, ratios_raw)
	    	print "Standardizing gene expression..."
	    	self.ratios_standardized = standardizeRatios(self, self.ratios)
	    	print "Inserting into row_info collection"
	    	self.row2id = get_row2id( self, self.ratios_standardized )
	    	self.row_info_collection = insert_row_info( self, self.row2id, row_annot )
	    	print "Inserting into col_info collection"
	    	self.col2id = get_col2id( self, self.ratios_standardized )
	    	self.col_info_collection = insert_col_info( self, self.col2id, col_annot )
	    	print "Inserting into ensemble_info collection"
	    	self.ensemble_info_collection = insert_ensemble_info( self, self.db_files, self.run2id, self.row2id, self.col2id )
	    	self.motif2gre = loadGREMap( self, gre2motif )

    	def checkRuns( self, db_files ):
		"""make sure the runs have data!!!"""
		to_keep = []
	    	for i in db_files:
	    		conn = sqlite3.connect( i )
		    	c = conn.cursor()
		    	c.execute("SELECT * FROM run_infos;")
		    	run_info = c.fetchone()
	    		if run_info is None:
	    			pass
	    		else:
	    			to_keep.append( i )
	    	return to_keep

    	def get_run2id( self, dbfiles ):
		"""make run2id"""
		runs = []
	    	for i in dbfiles:
	    		runs.append(i.split("/")[-2])
	    	row_info =  pd.DataFrame( zip( range( len( runs ) ), runs ), index = runs, columns = [ "run_id", "run_name"] )
	    	return row_info
	
	def loadRatios( self, file_in ):
		"""Loads ratios from individual cMonkey runs (unfinished) or single flat file (gzip compressed)."""
		if file_in == None:
			# compile from individual runs
			# do be done
			file_in = np.sort( np.array( glob.glob( e_dir + prefix + "???/ratios.tsv.gz" ) ) )
		else:
			print "Loading gene expression file from %s" % ratios_raw
			# load directly from gzip file
			ratios= pd.read_csv( gzip.open( file_in, 'rb' ), index_col=0, sep="\t" ) 
		return ratios

	def standardizeRatios( self, ratios ):
		"""compute standardized ratios (global). row standardized"""
		ratios_standardized = ratios.copy()
		zscore = lambda x: ( x - x.mean() ) / x.std()
		for row in ratios.iterrows():
			ratios_standardized.loc[ row[0] ] = zscore( row[1] )
		return ratios_standardized

	def get_row2id( self, ratios_standardized ):
		"""make row2id and id2row dicts for lookup"""
	    	rows = ratios_standardized.index.values
	    	rows.sort()
	    	row_info =  pd.DataFrame( zip( range( len( rows ) ), rows ), index = rows, columns = [ "row_id", "egrin2_row_name"] )
	    	return row_info

	def insert_row_info( self, row_info, row_annot, left_on = None, right_on = None ):
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
	    	# load additional row info
	    	if  row_annot != None:
	    		# assumes default microbes online schema
	    		row_annot = pd.read_csv( gzip.open( row_annot, 'rb' ), sep="\t" )	
	    	if left_on == None:
	    		left_on="egrin2_row_name"
	    	if right_on == None:
	    		right_on="sysName"
	    	# join with row_annot
	    	row_table = pd.merge( self.row_info, row_annot, left_on=left_on, right_on=right_on )
	    	# write to mongoDB collection 
	    	
	    	row_info_collection = db.row_info
	    	# TODO:
	    	# Check whether documents are already present in the collection before insertion
	    	# In case where they are present, update them
	    	row_info_collection.insert( row_table.to_dict( outtype='records' ) )

	    	return row_info_collection

	def get_col2id( self, ratios_standardized ):
		"""make cond2id and id2cond dicts for lookup"""
	    	cols = ratios_standardized.columns.values
	    	cols.sort()
	    	col_info =  pd.DataFrame( zip( range( len( cols ) ), cols ), index = cols, columns = [ "col_id", "egrin2_col_name"] )
	    	return col_info

    	def insert_col_info( self, col_info, col_annot ):
		"""
		Insert col_info into mongoDB database


	    	example queries
	    	------------------------------
	    	for i in col_info_collection.find( { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ):
	    		print i["egrin2_col_name"]
	    	for i in col_info_collection.find( { "$and": [ { "additional_info.name": "strain", "additional_info.value":  { "$regex": "MG1655" } }, { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ] } ):
	    		print i["egrin2_col_name"]
	    	for i in col_info_collection.find( { "$and": [ { "additional_info.value":   "0.0000303"  }, { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"  } ] } ):
	    		print i["egrin2_col_name"]

		"""
	    	# load additional col info
	    	
	    	if  col_annot != None:
	    		# assumes default microbes online schema
	    		col_annot = pd.read_csv( gzip.open( col_annot, 'rb' ), sep="\t" )	
	    	
	    	col_table = pd.merge( self.col_info, col_annot, left_on="egrin2_col_name", right_on="experiment_name" )
	    	col_info_4_mongoDB = []
	    	for i in range( 0, len( col_info.egrin2_col_name ) ):
	    		col_info_4_mongoDB.append( condInfo2Dict( 0, col_table, col_info.egrin2_col_name[i] ) )
	    	col_info_collection = db.col_info
	    	# TODO:
	    	# Check whether documents are already present in the collection before insertion
	    	# In case where they are present, update them
	    	col_info_collection.insert( col_info_4_mongoDB )

	    	return col_info_collection
	    	 	
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

	def assemble_ensemble_info( self, db_file, run2id, row2id, col2id ):
		"""Create python ensemble_info dictionary for bulk import into MongoDB collections"""
		run_name = db_file.split("/")[-2]
		conn = sqlite3.connect(db_file)
	    	c = conn.cursor()
	    	c.execute("SELECT start_time, finish_time, num_iterations, organism, species, num_rows, num_columns, num_clusters, git_sha FROM run_infos;")
	    	run_info = c.fetchone()
	    	c.execute("SELECT name FROM row_names;")
	    	rows = [ row2id.loc[ str(i[0]) ].row_id for i in c.fetchall() ]
	    	c.execute("SELECT name FROM column_names;")
	    	cols = [ col2id.loc[ str(i[0]) ].col_id for i in c.fetchall() ]
	    	d = {
	    	"run_id": run2id.loc[run_name].run_id,
	    	"run_name": run_name,
	    	"start_time": str( run_info[0] ),
	    	"finish_time": str( run_info[1] ),
	    	"num_iterations": int( run_info[2] ),
	    	"organism": str( run_info[3] ),
	    	"species": str( run_info[4] ),
	    	"num_rows": int( run_info[5] ),
	    	"rows": rows,
	    	"num_columns": int( run_info[6] ),
	    	"cols": cols,
	    	"num_clusters": int( run_info[7] ),
	    	"git_sha": str( run_info[8] ),
	    	"added_to_ensemble": datetime.datetime.utcnow()
	    	}
	    	conn.close()
	    	return d

	def insert_ensemble_info( self, db_files, run2id, row2id, col2id ):
		"""Compile and insert ensemble_info collection into MongoDB collection"""
	    	to_insert = [ assemble_ensemble_info( self, i, run2id, row2id, col2id ) for i in db_files ]
	    	ensemble_info_collection = db.ensemble_info
	    	# TODO:
	    	# Check whether documents are already present in the collection before insertion
	    	# In case where they are present, update them
	    	ensemble_info_collection.insert( to_insert )
	    	return ensemble_info_collection

 	def loadGREMap( self, gre2motif ):
 		count = 1
 		mots = {}
 		with open(gre2motif, 'r') as f: 
 			for line in f:
 				for motif in line.strip("\n").split( "\t" ):
 					elements = motif.split("_")
 					elements[1] = int(elements[1])
 					elements[2] = int(elements[2])
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
				count = count + 1
 		return mots


	def assemble_bicluster_info_single( self, db_file, cluster, run2id, row2id, col2id, motif2gre, row_info_collection ):
		"""Create python ensemble_info dictionary for bulk import into MongoDB collections"""
		run_name = db_file.split("/")[-2]
		conn = sqlite3.connect(db_file)
	    	c = conn.cursor()
	    	c.execute("SELECT max(iteration) FROM cluster_stats;")
	    	last_run = c.fetchone()[0] # i think there is an indexing problem in cMonkey python!! 
	    	w = (cluster,last_run)
	    	c.execute("SELECT residual FROM cluster_stats WHERE cluster = ? AND iteration = ?;", w )
	    	residual = c.fetchone()[0]
	    	c.execute("SELECT name FROM row_members JOIN row_names ON row_members.order_num = row_names.order_num WHERE row_members.cluster = ? AND row_members.iteration = ?;", w )
	    	rows = [ row2id.loc[ str(i[0]) ].row_id for i in c.fetchall() ]
	    	c.execute("SELECT name FROM column_members JOIN column_names ON column_members.order_num = column_names.order_num WHERE column_members.cluster = ? AND column_members.iteration = ?;", w )
	    	cols = [ col2id.loc[ str(i[0]) ].col_id for i in c.fetchall() ]
	    	c.execute("SELECT motif_num FROM motif_infos WHERE cluster = ? AND iteration = ?;", w )
	    	motif_nums = [ i[0] for i in c.fetchall() ]
	    	d = {
	    	"run_id": run2id.loc[run_name].run_id,
	    	"rows":, rows,
	    	"columns": cols,
	    	"residual": residual,
	    	"motif": [ get_motif_info_single(self, c, last_run, run_name, cluster, i, motif2gre) for i in motif_nums ]
	    	}
	    	conn.close()
	    	return d

   	def get_motif_info_single( self, cursor, iteration, run_name, cluster, motif_num, motif2gre, row_info_collection):
   		w = (cluster, iteration, motif_num)
   		cursor.execute("SELECT seqtype, evalue FROM motif_infos WHERE cluster = ? AND iteration = ? AND motif_num = ?;", w )
	    	motif_data = cursor.fetchone()
	    	cursor.execute("SELECT meme_motif_sites.rowid FROM meme_motif_sites JOIN motif_infos ON meme_motif_sites.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ?;", w )
	    	meme_site = [ i[0] for i in cursor.fetchall()]
   		d = {
	    	"gre_id": motif2gre[run_name][cluster][motif_num],
	    	"motif_num": motif_num,
	    	"seqtype": motif_data[0],
	    	"evalue": motif_data[1],
	    	"meme_motif_site": [],
	    	"pwm": {}
	    	}
	    	return d

	def get_meme_motif_site_single( self, cursor, iteration, run_name, cluster, motif_num, rowid, row_info_collection ):
		w = (cluster, iteration, motif_num, rowid)
		cursor.execute("SELECT seq_name, reverse, start, pvalue, flank_left, seq, flank_right FROM meme_motif_sites JOIN motif_infos ON meme_motif_sites.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ? AND meme_motif_sites.rowid = ?;", w )
		data = cursor.fetchone()

		d = {
		"row_id": row_info_collection.find( { "accession" : data[0] } )[0]["row_id"], #translate from accession to row_id, require microbes online table
		"reverse": data[1],
		"chr": row_info_collection.find( { "accession" : data[0] } )[0]["scaffoldId"],
		"start": data[3],
		"pvalue": data[4],
		"flank_left": data[5],
		"seq": data[6],
		"flank_right": data[7]
		}
		return(d)



