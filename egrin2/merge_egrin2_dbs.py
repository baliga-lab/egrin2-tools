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


# ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz"
# col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz"
# db_file = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/eco-out-001/cmonkey_run.db"

# delete egrin2 db. run on console

# sql2mongoDB( ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz",  col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz", ncbi_code = "511145")

class sql2mongoDB:
    
	def __init__( self, e_dir = None, prefix = None, gene_info = None, ratios_raw = None, gre2motif = None, col_annot = None, ncbi_code = None, dbname = None ):
		
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
		self.db = client[dbname]

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
			self.gre2motif = e_dir + "out.mot_metaclustering.txt.I45.txt"
		else:
			self.gre2motif = gre2motif

	    	self.db_files = np.sort( np.array( glob.glob( e_dir + prefix + "???/cmonkey_run.db" ) ) ) # get all cmonkey_run.db files
	    	
	    	if ncbi_code == None:
			# try to find ncbi_code in cmonkey_run.db
			conn = sqlite3.connect( self.db_files[0] )
		    	c = conn.cursor()
		    	try:
		    		c.execute("SELECT ncbi_code FROM run_infos;")
		    		self.ncbi_code = c.fetchone()[0]
		    	except sqlite3.Error as e:
		    		print "Could not find NCBI Genome ID in cmonkey_run.db:", e.args[0]
    		else:
    			self.ncbi_code = ncbi_code

    		self.ratios_raw = ratios_raw
    		self.col_annot = col_annot

    		return None

	def dlfile( self, url, save_name=None):
		# Open the url
		if save_name==None:
			save_name = "tmp"
		try:
			f = urlopen(url)
			print "downloading " + url

			# Open our local file for writing
			with open(os.path.basename(save_name), "wb") as local_file:
			    local_file.write(f.read())

		#handle errors
		except HTTPError, e:
			print "HTTP Error:", e.code, url
		except URLError, e:
			print "URL Error:", e.reason, url
    	
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

	def loadGenome (self, ncbi_code):
		# download genome from microbes online. store in MongoDB collection
		# 
		url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId="+ncbi_code+";export=genome"
		save_name=ncbi_code+"_genome.fa"
		self.dlfile( url, save_name)
		
		with open( save_name, 'r' ) as f:
			fasta_sequences = SeqIO.parse( f ,'fasta')
			for fasta in fasta_sequences:
	   			d = {
				"scaffoldId": fasta.id,
				"NCBI_RefSeq": fasta.description.split(" ")[1],
				"NCBI_taxonomyId": fasta.description.split(" ")[-1],
				"sequence": str(fasta.seq)
				}
		genome_collection = self.db.genome
	    	# TODO:
	    	# Check whether documents are already present in the collection before insertion
	    	# In case where they are present, update them
	    	genome_collection.insert( d )

	    	return genome_collection
	
	def loadRatios( self, file_in ):
		"""Loads ratios from individual cMonkey runs (unfinished) or single flat file (gzip compressed)."""
		if file_in == None:
			# compile from individual runs
			# do be done
			file_in = np.sort( np.array( glob.glob( e_dir + prefix + "???/ratios.tsv.gz" ) ) )
		else:
			print "Loading gene expression file from %s" % file_in
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

	def insert_row_info( self, ncbi_code, row_info, left_on = None, right_on = None ):
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
		if ncbi_code != None:
			print "Downloading gene information for NCBI taxonomy ID:", ncbi_code
		    	url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId="+ncbi_code+";export=tab"
			save_name=ncbi_code+"_geneInfo.tab"
			self.dlfile(url,save_name)
			
			with open( save_name, 'r' ) as f:
				row_annot = pd.read_csv( f, sep="\t" )
		    	if left_on == None:
		    		left_on="egrin2_row_name"
		    	if right_on == None:
		    		right_on="sysName"
		    	# join with row_annot
		    	row_table = pd.merge( row_info, row_annot, left_on=left_on, right_on=right_on )
		    	# write to mongoDB collection 
		    	
		    	row_info_collection = self.db.row_info
		    	# TODO:
		    	# Check whether documents are already present in the collection before insertion
		    	# In case where they are present, update them
		    	row_info_collection.insert( row_table.to_dict( outtype='records' ) )
	    	else:
	    		print "WARNING: could not fetch additional gene information for NCBI taxonomy ID:", ncbi_code
	    		row_info_collection = self.db.row_info
		    	# TODO:
		    	# Check whether documents are already present in the collection before insertion
		    	# In case where they are present, update them
		    	row_info_collection.insert( row_info.to_dict( outtype='records' ) )

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
			print i

		"""
		# load additional col info

		if  col_annot != None:
			# assumes default microbes online schema
			col_annot = pd.read_csv( gzip.open( col_annot, 'rb' ), sep="\t" )	

		col_table = pd.merge( col_info, col_annot, left_on="egrin2_col_name", right_on="experiment_name" )
		col_info_4_mongoDB = []
		for i in range( 0, len( col_info.egrin2_col_name ) ):
			col_info_4_mongoDB.append( self.condInfo2Dict( col_table, col_info.egrin2_col_name[i] ) )
		col_info_collection = self.db.col_info
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
	    	to_insert = [ self.assemble_ensemble_info( i, run2id, row2id, col2id ) for i in db_files ]
	    	ensemble_info_collection = self.db.ensemble_info
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

	def insert_bicluster_info( self, db_file, run2id, row2id, col2id, motif2gre, row_info_collection ): 
		"""Find all biclusters in a cMonkey run, process and add as documents to bicluster collection

		example queries
		------------------------------

		bicluster_info_collection.find({"rows":{"$all":[26,27]}}).count()
		"""
		# Get all biclusters from cmonkey run
		conn = sqlite3.connect(db_file)
	    	c = conn.cursor()
	    	c.execute("SELECT max(iteration) FROM cluster_stats;")
	    	last_run = c.fetchone()[0] # i think there is an indexing problem in cMonkey python!! 
	    	w = (last_run,)
	    	c.execute("SELECT cluster FROM cluster_stats WHERE iteration = ?;",w)
		biclusters = [self.assemble_bicluster_info_single( db_file, c, last_run, i[0], run2id, row2id, col2id, motif2gre, row_info_collection ) for i in c.fetchall()]
		bicluster_info_collection = self.db.bicluster_info
	    	# TODO:
	    	# Check whether documents are already present in the collection before insertion
	    	# In case where they are present, update them
	    	bicluster_info_collection.insert( biclusters )
	    	return None

	def assemble_bicluster_info_single( self, db_file, cursor, iteration, cluster, run2id, row2id, col2id, motif2gre, row_info_collection ):
		"""Create python ensemble_info dictionary for bulk import into MongoDB collections"""
		#print cluster
		run_name = db_file.split("/")[-2]
	    	w = (cluster,iteration)
	    	cursor.execute("SELECT residual FROM cluster_stats WHERE cluster = ? AND iteration = ?;", w )
	    	residual = cursor.fetchone()[0]
	    	cursor.execute("SELECT name FROM row_members JOIN row_names ON row_members.order_num = row_names.order_num WHERE row_members.cluster = ? AND row_members.iteration = ?;", w )
	    	rows = [ row2id.loc[ str(i[0]) ].row_id for i in cursor.fetchall() ]
	    	cursor.execute("SELECT name FROM column_members JOIN column_names ON column_members.order_num = column_names.order_num WHERE column_members.cluster = ? AND column_members.iteration = ?;", w )
	    	cols = [ col2id.loc[ str(i[0]) ].col_id for i in cursor.fetchall() ]
	    	cursor.execute("SELECT motif_num FROM motif_infos WHERE cluster = ? AND iteration = ?;", w )
	    	motif_nums = [ i[0] for i in cursor.fetchall() ]
	    	d = {
	    	"run_id": run2id.loc[run_name].run_id,
	    	"rows": rows,
	    	"columns": cols,
	    	"residual": residual,
	    	"motif": [ self.get_motif_info_single(cursor, iteration, run_name, cluster, i, motif2gre, row_info_collection) for i in motif_nums ]
	    	}
	    	return d

	def get_motif_info_single( self, cursor, iteration, run_name, cluster, motif_num, motif2gre, row_info_collection):
		w = (cluster, iteration, motif_num)
		cursor.execute("SELECT seqtype, evalue FROM motif_infos WHERE cluster = ? AND iteration = ? AND motif_num = ?;", w )
		motif_data = cursor.fetchone()
		cursor.execute("SELECT meme_motif_sites.rowid FROM meme_motif_sites JOIN motif_infos ON meme_motif_sites.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ?;", w )
		meme_site = [ i[0] for i in cursor.fetchall()]
		cursor.execute("SELECT row FROM motif_pssm_rows JOIN motif_infos ON motif_pssm_rows.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ?;", w )
		rows = [ i[0] for i in cursor.fetchall()]

		try:
			gre_id = motif2gre[run_name][cluster][motif_num] 
		except:
			gre_id = "NaN"

		d = {
		"gre_id": gre_id,
		"motif_num": motif_num,
		"seqtype": motif_data[0],
		"evalue": motif_data[1],
		"meme_motif_site": [self.get_meme_motif_site_single( cursor, iteration, run_name, cluster, motif_num, i, row_info_collection ) for i in meme_site],
		"pwm": [self.get_pwm_single( cursor, iteration, run_name, cluster, motif_num, i, row_info_collection ) for i in rows]
		}
		return d

	def get_meme_motif_site_single( self, cursor, iteration, run_name, cluster, motif_num, rowid, row_info_collection ):
		w = (str(rowid),)
		cursor.execute("SELECT seq_name, reverse, start, pvalue, flank_left, seq, flank_right FROM meme_motif_sites WHERE rowid = ?", w )
		data = cursor.fetchone()

		# try to match accession to row_id
		try:
			row_id = row_info_collection.find( { "accession" : data[0] } )[0]["row_id"] #translate from accession to row_id, requires microbes online
		except:
			row_id = "NaN"

		try:
			scaffoldId = row_info_collection.find( { "accession" : data[0] } )[0]["scaffoldId"]
		except:
			scaffoldId = "NaN"

		d = {
		"row_id": row_id,
		"reverse": data[1],
		"scaffoldId": scaffoldId,
		"start": data[2],
		"pvalue": data[3],
		"flank_left": data[4],
		"seq": data[5],
		"flank_right": data[6]
		}
		return d

	def get_pwm_single( self, cursor, iteration, run_name, cluster, motif_num, row, row_info_collection ):
		w = (cluster, iteration, motif_num, row)
		cursor.execute("SELECT row, a, c, g, t FROM motif_pssm_rows JOIN motif_infos ON motif_pssm_rows.motif_info_id = motif_infos.rowid WHERE motif_infos.cluster = ? AND motif_infos.iteration = ? AND motif_infos.motif_num = ? AND motif_pssm_rows.row = ?;", w )
		data = cursor.fetchone()

		d = {
		"row": data[0],
		"a": data[1],
		"c": data[2],
		"g": data[3],
		"t": data[4]
		}
		return d

	def compile( self ):
		# print "Compiling EGRIN2 ensemble..."  
		self.db_files = self.checkRuns( self.db_files )
	    	self.run2id = self.get_run2id( self.db_files )
	    	print "Downloading genome information for NCBI taxonomy ID:", self.ncbi_code
		self.genome_collection = self.loadGenome( self.ncbi_code )
	    	self.ratios = self.loadRatios( self.ratios_raw)
	    	print "Standardizing gene expression..."
	    	self.ratios_standardized = self.standardizeRatios( self.ratios)
	    	print "Inserting into row_info collection"
	    	self.row2id = self.get_row2id( self.ratios_standardized )
	    	self.row_info_collection = self.insert_row_info( self.ncbi_code, self.row2id )
	    	print "Inserting into col_info collection"
	    	self.col2id = self.get_col2id( self.ratios_standardized )
	    	self.col_info_collection = self.insert_col_info( self.col2id, self.col_annot )
	    	print "Inserting into ensemble_info collection"
	    	self.ensemble_info_collection = self.insert_ensemble_info( self.db_files, self.run2id, self.row2id, self.col2id )
	    	self.motif2gre = self.loadGREMap( self.gre2motif )
	    	print "Inserting into bicluster collection"
	    	for i in self.db_files:
	    		print i
	    		self.insert_bicluster_info( i, self.run2id, self.row2id, self.col2id, self.motif2gre, self.row_info_collection )
	    	return None



