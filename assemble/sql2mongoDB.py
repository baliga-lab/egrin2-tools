#!/usr/bin/env python

"""
Initialize  mongo database from individual cMonkey 
runs (SQLite) plus some additional tables and run_id column, 
including motifs and motif clusters
"""

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
from bson.objectid import ObjectId
import gridfs
from Bio import SeqIO


# ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz"
# col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz"
# db_file = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/eco-out-001/cmonkey_run.db"

# delete egrin2 db. run on console

# sql2mongoDB( ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz",  col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz", ncbi_code = "511145")

class sql2mongoDB:
    
	def __init__( self, organism = None, host = None, port = None, e_dir = None, targetdir = None, prefix = None,ratios_raw = None, gre2motif = None, col_annot = None, ncbi_code = None, dbname = None , db_run_override = None, genome_file = None, row_annot = None, row_annot_match_col = None ):
		
		# connect to database
		# make sure mongodb is running
		# retvalue = os.system("nohup mongod --port 27017 --quiet &")

		if host is None:
			host = "localhost"

		if port is None:
			port = 27017

		if organism is None:
			print "Requires an organism code, e.g. eco for E. coli"
			return None
		else:
			self.organism = organism

		try:
			client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' ) 
			print "Connected to MongoDB"
		except pymongo.errors.ConnectionFailure, e:
			print "Could not connect to MongoDB: %s" % e

		if dbname is None:
			self.dbname = self.organism + "_db"
		else:
			self.dbname = dbname

		if self.dbname in client.database_names():
			print "WARNING: %s database already exists!!!" % self.dbname
		else:
			print "Initializing MongoDB database: %s" % self.dbname
		self.db = client[self.dbname]

		# get db files in directory
		if prefix is None:
			self.prefix = organism+'-out-'
		else:
			self.prefix = prefix
		if e_dir is None:
			self.e_dir = './'
		else:
			self.e_dir = e_dir
		if targetdir is None:
			self.targetdir = './'
		else:
			self.targetdir = targetdir
		if gre2motif == None:
			# default file name?
			self.gre2motif = self.e_dir + "out.mot_metaclustering.txt.I45.txt"
		else:
			self.gre2motif = gre2motif

	    	self.db_files = np.sort( np.array( glob.glob( self.e_dir + self.prefix + "???/cmonkey_run.db" ) ) ) # get all cmonkey_run.db files
	    	self.db_run_override = db_run_override
	    	
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
    		self.genome_file = genome_file
    		self.row_annot = row_annot
    		self.row_annot_match_col = row_annot_match_col

    		if len(self.db_files) < 1:
	    		print "I cannot find any cMonkey SQLite databases in the current directory: %s\nMake sure 'e_dir' variable points to the location of your cMonkey-2 ensemble results." % os.getcwd()

    		return None

	def dlfile( self, url, save_name=None):
		# Open the url
		if save_name==None:
			save_name = "tmp"
		count = 1
		while count < 5:
			try:
				f = urlopen(url)
				print "downloading " + url

				# Open our local file for writing
				with open(os.path.basename(save_name), "wb") as local_file:
				    local_file.write(f.read())
				break
			#handle errors
			except HTTPError, e:
				print "HTTP Error:", e.code, url
				count = count + 1
				print "Trying to connect again. Attempt %s of 5" % count
			except URLError, e:
				print "URL Error:", e.reason, url
				print "Trying to connect again. Attempt %s of 5" % count
				count = count + 1
    	
	def checkRuns( self, db_files, db_run_override, db ):
		"""make sure the runs have data!!!"""
		to_keep = []
		ensemble_info_collection = db.ensemble_info
		for i in db_files:
			conn = sqlite3.connect( i )
			c = conn.cursor()
			c.execute("SELECT * FROM run_infos;")
			run_info = c.fetchone()
			if run_info is None:
				pass
			else:
				if db_run_override == None:
					# do not include runs that are already in the database
					# check for existence
					run_name = i.split("/")[-2]
					if ensemble_info_collection.find( { "run_name": run_name } ).count() > 0:
						pass
					else:
						to_keep.append( i )
				else:
					to_keep.append( i )
		return to_keep

	def get_run2id( self, dbfiles, db ):
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
					run_id.append(max(run_id)+1)
				else:
					run_id.append(0)
		run_info =  pd.DataFrame( zip( run_id, run_name ), index = run_name, columns = [ "run_id", "run_name"] )
		return run_info

	def check4existence( self, collection, document, key1 = None, value1 = None, key2 = None, value2 = None ):
		if key1 == None:
			d_check = collection.find( document ).count()
		else:
			if key2 == None:
				d_check = collection.find( { key1 : value1 } ).count()
			else:
				d_check = collection.find( { key1 : value1, key2 : value2 } ).count()
		if d_check == 0:
			return document

	def loadGenome (self, ncbi_code, genome_file):

		if genome_file == None:
			print "No custom genome annotation file supplied by 'genome_file' parameter. Attempting automated download from MicrobesOnline"
			# download genome from microbes online. store in MongoDB collection
			# 
			url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId="+ncbi_code+";export=genome"
			save_name=ncbi_code+"_genome.fa"
			self.dlfile( url, save_name)
			
			seqs_b = []
			with open( save_name, 'r' ) as f:
				fasta_sequences = SeqIO.parse( f ,'fasta')
				for fasta in fasta_sequences:
		   			seqs_b.append( {
					"scaffoldId": fasta.id,
					"NCBI_RefSeq": fasta.description.split(" ")[1],
					"NCBI_taxonomyId": fasta.description.split(" ")[-1],
					"sequence": str(fasta.seq)
					} )
		else:
			# TBD
			seqs_b = []
		
		genome_collection = self.db.genome
	    
	    	# Check whether documents are already present in the collection before insertion
	    	seqs_f = filter( None, [ self.check4existence( genome_collection, i ) for i in seqs_b ] )

	    	print "%s new records to write" % len( seqs_f )
	    	if len(seqs_f) > 0:
	    		genome_collection.insert( seqs_f )

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

	def get_row2id( self, ratios_standardized, db ):
		"""make row2id and id2row dicts for lookup"""
		row_info_collection = db.row_info
		row_name = []
		row_id = []
		for i in row_info_collection.find():
			row_name.append( i["egrin2_row_name"] )
			row_id.append( i["row_id"] )
		for i in ratios_standardized.index.values:
			if i not in row_name:
				row_name.append( i )
			if len(row_id) > 0:
				row_id.append( max(row_id)+1 )
			else:
				row_id.append(0)
		row_info =  pd.DataFrame( zip(row_id, row_name), index = row_name, columns = [ "row_id", "egrin2_row_name"] )
		return row_info

	def insert_row_info( self, ncbi_code, row_info, row_annot, row_annot_match_col ):
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
			print "Row annotation file not suppled as 'row_annot' parameter. Attempting automated download from MicrobesOnline"
			if ncbi_code != None:
				print "Downloading gene information for NCBI taxonomy ID:", ncbi_code
			    	url = "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId="+ncbi_code+";export=tab"
				save_name=ncbi_code+"_geneInfo.tab"
				self.dlfile(url,save_name)
				
				with open( save_name, 'r' ) as f:
					row_annot = pd.read_csv( f, sep="\t" )
			    	left_on="egrin2_row_name"
			    	if row_annot_match_col == None:
			    		row_annot_match_col="sysName"
			    	# join with row_annot
			    	row_table = pd.merge( row_info, row_annot, left_on=left_on, right_on=row_annot_match_col )	
		    	else:
		    		print "WARNING: could not fetch additional gene information for NCBI taxonomy ID:", ncbi_code
			    	# TODO:
			    	# Check whether documents are already present in the collection before insertion
			    	# In case where they are present, update them
			    	row_table = row_info
	    	else:
			row_annot = pd.read_csv( open( row_annot, 'rb' ), sep="\t" )	
			# join with row_annot
		    	row_table = pd.merge( row_info, row_annot, left_on=left_on, right_on=row_annot_match_col )	

		# write to mongoDB collection 
		row_info_collection = self.db.row_info
		# Check whether documents are already present in the collection before insertion
		d = row_table.to_dict( orient='records' )
	    	d_f = filter( None, [ self.check4existence( row_info_collection, i ) for i in d ] )

	    	print "%s new records to write" % len( d_f )
	    	
	    	if len(d_f) > 0:
	    		row_info_collection.insert( d_f )

	    	return row_info_collection

	def get_col2id( self, ratios_standardized, db ):
		"""make cond2id and id2cond dicts for lookup"""
		col_info_collection = db.col_info
		col_name = []
		col_id = []
		for i in col_info_collection.find():
			col_name.append(i["egrin2_col_name"])
			col_id.append(i["col_id"])
		for i in ratios_standardized.columns.values:
	    		if i not in col_name:
	    			col_name.append( i )
	    			if len(col_id) > 0:
					col_id.append(max(col_id)+1)
				else:
					col_id.append(0)
	    	col_info =  pd.DataFrame( zip( col_id, col_name ), index = col_name, columns = [ "col_id", "egrin2_col_name"] )
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
		
		# write to mongoDB collection 
		col_info_collection = self.db.col_info
		
		# Check whether documents are already present in the collection before insertion
	    	d_f = filter( None, [ self.check4existence( col_info_collection, i ) for i in col_info_4_mongoDB ] )

	    	print "%s new records to write" % len( d_f )
	    	
	    	if len(d_f) > 0:
	    		col_info_collection.insert( d_f )

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

	def insert_gene_expression( self, db, row2id, col2id, ratios, ratios_standardized ):
		"""
		Insert gene_expression into mongoDB database


		example queries
		------------------------------

		"""
		exp_data = []
		counter = 0
		for i in ratios.index.values:
			if counter%200 == 0:
				print "%s percent done" % round( ( float(counter)/ratios.shape[0] )*100,1 )
			for j in ratios.columns.values:
				exp_data.append(
					{
				    	"row_id": row2id.loc[i].row_id,
				    	"col_id": col2id.loc[j].col_id,
			 		"normalized_expression": ratios.loc[i,j],
			 		"standardized_expression": ratios_standardized.loc[i,j]
				    	} )
			counter = counter + 1

		# write to mongoDB collection 
		gene_expression_collection = db.gene_expression
		
		# Check whether documents are already present in the collection before insertion
	    	d_f = filter( None, [ self.check4existence( gene_expression_collection, i ) for i in exp_data ] )

	    	print "%s new records to write" % len( d_f )
	    	
	    	if len(d_f) > 0:
	    		gene_expression_collection.insert( d_f )

		return gene_expression_collection

	def assemble_ensemble_info( self, db_file, run2id, row2id, col2id ):
		"""Create python ensemble_info dictionary for bulk import into MongoDB collections"""  
		run_name = db_file.split("/")[-2]
		print "Assembling run info for cMonkey run: %s" % run_name
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

	def insert_ensemble_info( self, db_files, db, run2id, row2id, col2id ):
		"""Compile and insert ensemble_info collection into MongoDB collection"""
	    	to_insert = [ self.assemble_ensemble_info( i, run2id, row2id, col2id ) for i in db_files ]
	    	ensemble_info_collection = db.ensemble_info
	    	
	    	# Check whether documents are already present in the collection before insertion
	    	d_f = filter( None, [ self.check4existence( ensemble_info_collection, i, "run_name", i["run_name"] ) for i in to_insert ] )

	    	print "%s new records to write" % len( d_f )
	    	
	    	if len(d_f) > 0:
	    		ensemble_info_collection.insert( d_f )

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

	def insert_bicluster_info( self, db, e_dir, db_file, run2id, row2id, col2id, motif2gre, row_info_collection ): 
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
		biclusters = [self.assemble_bicluster_info_single( db, e_dir, db_file, c, last_run, i[0], run2id, row2id, col2id, motif2gre, row_info_collection ) for i in c.fetchall()]
		bicluster_info_collection = self.db.bicluster_info
	    	# Check whether documents are already present in the collection before insertion
	    	d_f = filter( None, [ self.check4existence( bicluster_info_collection, i, "run_id", i["run_id"], "cluster", i["cluster"] ) for i in biclusters ] )

	    	print "%s new records to write" % len( d_f )
	    	
	    	if len(d_f) > 0:
	    		bicluster_info_collection.insert( d_f )
	    	
	    	return bicluster_info_collection

	def assemble_bicluster_info_single( self, db, e_dir, db_file, cursor, iteration, cluster, run2id, row2id, col2id, motif2gre, row_info_collection ):
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
	    	"cluster": cluster,
	    	"rows": rows,
	    	"columns": cols,
	    	"residual": residual,
	    	"motif": [ self.get_motif_info_single( db, e_dir, cursor, iteration, run_name, cluster, i, motif2gre, row_info_collection) for i in motif_nums ]
	    	}
	    	return d

	def get_motif_info_single( self, db, e_dir, cursor, iteration, run_name, cluster, motif_num, motif2gre, row_info_collection):
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
		# do not store this info
		# anb 01/22/2015
		# "flank_left": data[4],
		# "seq": data[5],
		# "flank_right": data[6]
		"pvalue": data[3],
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

	def tmp_fcn( self, i, e_dir ):
		print i.cluster, i.run_id, e_dir
		return None

	def assemble_fimo( self ):

		def get_fimo_scans_single( i, db, e_dir, run2id ):
			cluster = i.cluster
			run_name =run2id[ run2id[ "run_id" ]==i.run_id ][ "run_name" ][ 0 ]
			cluster_id = i._id
			try:
				# get all fimo scans in the dir
				f = e_dir + run_name + "/fimo-outs/fimo-out-" + "%04d" % (cluster,) + ".bz2" 
				#print motif_num

				fimo = pd.read_csv( f, sep="\t", compression = "bz2" )
				# change sequence_name to scaffoldId
				fimo.rename(columns={'matched sequence': 'matched_sequence', 'sequence name': 'scaffoldId', "#pattern name": "motif_num"}, inplace=True)

				trans_d = {}
				for i in np.unique(fimo.scaffoldId):
					NCBI_RefSeq = "_".join(i.split('.')[-2].split("_")[::-1][0:2][::-1])
					scaffoldId = db.genome.find_one( { "NCBI_RefSeq": NCBI_RefSeq } )["scaffoldId"]
					trans_d[i] = scaffoldId

				trans_v = [trans_d[i] for i in fimo.scaffoldId.values]
				fimo.scaffoldId = trans_v
				fimo["cluster_id"] = cluster_id
				#fimo["cluster"] = cluster

				# only keep specific columns
				fimo = fimo.loc[ : , [ 'scaffoldId', 'start', 'stop', 'strand', 'score', 'p-value', 'in_coding_rgn', 'cluster_id' ] ]

				d_f = fimo.to_dict( orient='records' )

				db.fimo.insert( d_f )

		    		return None
			except Exception:
				return None

		# get all biclusters
		bcs = pd.DataFrame( list( self.db.bicluster_info.find( {}, { "cluster" : 1, "run_id": 1 } ) ) )
		tmp = bcs.apply( get_fimo_scans_single, axis=1, db = self.db, e_dir = self.e_dir, run2id = self.run2id  )
		return None

	def mongoDump( self, db, outfile ):
		"""Write contents from MongoDB instance to binary file"""
		outfile = os.path.abspath( os.path.join( self.targetdir,outfile ) )
		sys_command = "mongodump --db " + db + " --out " + outfile + "--host " + self.host + "--port " +self.port 
		os.system(sys_command)

	def mongoRestore( self, db, infile ):
		"""Read contents of binary MongoDB dump into MongoDB instance"""
		sys_command = "mongorestore --db " + db + " " + infile + "--host " + self.host + "--port " +self.port 
		os.system(sys_command)

	def compile( self ):
		"""Compile EGRIN2 ensemble"""
		# print "Compiling EGRIN2 ensemble..."  
		self.db_files = self.checkRuns( self.db_files, self.db_run_override, self.db )
	    	self.run2id = self.get_run2id( self.db_files, self.db )
	    	
	    	print "Downloading genome information for NCBI taxonomy ID:", self.ncbi_code
		self.genome_collection = self.loadGenome( self.ncbi_code, self.genome_file )
	    	self.expression = self.loadRatios( self.ratios_raw)
	    	
	    	print "Standardizing gene expression..."
	    	self.expression_standardized = self.standardizeRatios( self.expression)
	    	
	    	print "Inserting into row_info collection"
	    	self.row2id = self.get_row2id( self.expression_standardized, self.db )
	    	self.row_info_collection = self.insert_row_info( self.ncbi_code, self.row2id, self.row_annot, self.row_annot_match_col )
	    	
	    	print "Inserting into col_info collection"
	    	self.col2id = self.get_col2id( self.expression_standardized, self.db )
	    	self.col_info_collection = self.insert_col_info( self.col2id, self.col_annot )
	    	
	    	print "Inserting gene expression into database"
	    	self.gene_expression_collection = self.insert_gene_expression( self.db, self.row2id, self.col2id, self.expression, self.expression_standardized )
	    	
	    	print "Inserting into ensemble_info collection"
	    	self.ensemble_info_collection = self.insert_ensemble_info( self.db_files, self.db, self.run2id, self.row2id, self.col2id )
	    	self.motif2gre = self.loadGREMap( self.gre2motif )
	    	
	    	print "Inserting into bicluster collection"
	    	for i in self.db_files:
	    		print i
	    		self.bicluster_info_collection = self.insert_bicluster_info( self.db, self.e_dir, i, self.run2id, self.row2id, self.col2id, self.motif2gre, self.row_info_collection )
		
		print "Indexing bicluster collection"
		self.bicluster_info_collection.ensure_index( "rows" )
		self.bicluster_info_collection.ensure_index( "columns" )
		self.bicluster_info_collection.ensure_index( "motif.gre_id" ) 
		
		print "Inserting into fimo collection. This might take awhile..."
		self.assemble_fimo( )
		
		print "Indexing fimo collection"
		self.db.fimo.ensure_index( [ ( "scaffoldId", pymongo.ASCENDING ), ( "start", pymongo.ASCENDING ), ( "stop", pymongo.ASCENDING ),( "p-value", pymongo.ASCENDING ), ( "cluster_id", pymongo.ASCENDING ) ] ) 


    		outfile =  self.prefix + str(datetime.datetime.utcnow()).split(" ")[0] + ".mongodump"
		
		print "Writing EGRIN2 MongoDB to %s" % self.targetdir + outfile   		
    		self.mongoDump( self.dbname, outfile )
	    	return None

if __name__ == '__main__':
	self = sql2mongoDB( ratios_raw = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz",  col_annot = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz", ncbi_code = "511145", e_dir = "/Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/", organism="eco", host = "primordial")


