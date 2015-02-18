#!/usr/bin/env python

"""Tools for querying EGRIN2.0 MongoDB."""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks", "David Reiss"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import random

# $ hg clone ssh://hg@bitbucket.org/djcbeach/monary ./monary
# $ cd ./monary && python setup.py install
#from monary import Monary
# can't get monary to connect to a remote host!!

from pymongo import MongoClient
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import itertools
from bson.code import Code
import matplotlib.pyplot as plt
import plotly.plotly as py
from plotly.graph_objs import *
import colorbrewer as cb
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

from assemble.resample import *

def rsd( vals ):
	return abs( np.std( vals ) / np.mean( vals ) )

def check_colResamples( col, n_rows, n_resamples, host="localhost", port=27017, db="" ):
	# connect to db
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' ) 
	if client[db].col_resample.find_one( { "n_rows": n_rows, "col_id": col, "resamples": { "$gte": n_resamples } } ) is None:
		client.close()
		return col
	client.close()

def findMatch( x, df, return_field ):
	print x
	"""Find which 'df' element x matches. Return appropriate translation"""
	# find matching row
	counter = 0
	while ( x in df.iloc[ counter ].tolist() ) == False:
		counter = counter + 1
	return df.iloc[ counter ][ return_field ]

def row2id( row, host="localhost", port=27017, db="",  verbose = False, return_field = "row_id" ):
	"""Check name format of rows. If necessary, translate."""
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = list( client[db].row_info.find( { "$or": [ { "row_id": row }, { "egrin2_row_name": row }, { "GI": row }, { "accession": row }, { "name": row }, { "sysName": row } ] } ) )
	client.close()
	if len( query ) == 1: 
		try:
			row = query[ 0 ][ return_field ]
		except Exception:
			row = query[ 0 ][ "row_id" ]
		return row
	elif len( query ) > 0:
		print "ERROR: Multiple rows match the row name: %s" % row
		if verbose:
			print query
		return None
	else:
		print "ERROR: Cannot identify row name: %s" % row
		return None

def row2id_batch( rows, host="localhost", port=27017, db="",  verbose = True, return_field = "row_id", input_type = None ):
	"""Check name format of rows. If necessary, translate."""

	if return_field == input_type:
		return rows
	
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = pd.DataFrame( list( client[db].row_info.find( { "$or": [ { "row_id": { "$in": rows } }, { "egrin2_row_name": { "$in": rows } }, { "GI": { "$in": rows } }, { "accession": {"$in": rows } }, { "name": { "$in": rows } }, { "sysName": { "$in": rows }  } ] }, { "row_id": 1, "egrin2_row_name": 1, "GI": 1, "accession" : 1, "name" : 1 } ) ) )

	if input_type in [  "row_id", "egrin2_row_name", "GI", "accession", "name", "sysName" ]:
		query = query.set_index(input_type)
		to_r = query.loc[ rows ][ return_field ].tolist()
	else:
		#try to match input_type automatically
		if verbose:
			print "Reverting to translation by single matches. Defining 'input_type' will dramatically speed up query."
		to_r= [ row2id( x, host, port, db, return_field ) for x in rows ]

	client.close()
	
	return to_r
	

def col2id( col, host="localhost", port=27017, db="",  verbose = False, return_field = "col_id" ):
	"""Check name format of rows. If necessary, translate."""

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = list( client[db].col_info.find( { "$or": [ { "col_id": col }, { "egrin2_col_name": col } ] } ) )
	client.close()
	if len( query ) == 1: 
		try:
			col = query[ 0 ][ return_field ]
		except Exception:
			col = query[ 0 ][ "col_id" ]
		return col
	elif len( query ) > 0:
		print "ERROR: Multiple cols match the col name: %s" % col
		if verbose:
			print query
		return None
	else:
		print "ERROR: Cannot identify col name: %s" % col
		return None


def col2id_batch( cols, host="localhost", port=27017, db="",  verbose = True, return_field = "col_id", input_type = None ):
	"""Check name format of rows. If necessary, translate."""

	if return_field == input_type:
		return cols

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = pd.DataFrame( list( client[db].col_info.find( { "$or": [ { "col_id": { "$in": cols } }, { "egrin2_col_name": { "$in":cols } } ] }, { "col_id": 1, "egrin2_col_name":1 } ) ) )

	if input_type in [  "col_id", "egrin2_col_name" ]:
		query = query.set_index( input_type )
		to_r = query.loc[ cols ][ return_field ].tolist()
	else:
		#try to match input_type automatically
		if verbose:
			print "Reverting to translation by single matches. Defining 'input_type' will dramatically speed up query."
		to_r= [ col2id( x, host, port, db, return_field ) for x in cols ]
	
	client.close()

	return to_r

def colResamplePval( rows = None, row_type = None, cols = None, col_type = None, n_resamples = None, host = "localhost", port = 27017, db = "", standardized = None, sig_cutoff = None, sort = True, add_override = False, n_jobs = 4, keepP = 0.1, verbose = True, col_outtype = "col_id" ):

	def empirical_pval( i, random_rsd, resamples ):
		for x in range( 0, len( i ) ):
			val = ( float ( sum( [ i.values[ 0 ] >= y for y in random_rsd[ i.index[0] ] ] ) ) / len( random_rsd[ i.index[0] ] ) ) * ( float( len( random_rsd[ i.index[0] ] ) ) / resamples[ i.index[0] ] )
			if val >= ( float( len( random_rsd[ i.index[0] ] ) ) / resamples[ i.index[0] ] ):
				#return ">= %f" % ( round( float( len( random_rsd[ "lowest_normalized" ] ) ) / random_rsd[ "resamples" ], 2 ) )
				return round( float( len( random_rsd[ i.index[0] ] ) ) / resamples[ i.index[0] ], 2 )
			elif val == 0:
				#return "< %f" % ( 1 / float( random_rsd[ "resamples" ] ) )
				return ( 1 / float( resamples[ i.index[0] ] ) )
			else:
				return val

	rows_o = rows
	rows = row2id_batch( rows, host, port, db, input_type = row_type )
	if len( rows ) == 0:
		print "Please provide an appropriately named array of rows"
		return None 
	
	cols_o = cols
	cols = col2id_batch( cols, host, port, db, input_type = col_type )
	if len( cols ) == 0:
		print "Please provide an appropriately named array of cols"
		return None

	if standardized is None:
		# compute RSD on standardized gene expression by default
		# other option 'False' for normalized (not standardized expression)
		standardized = True

	if sig_cutoff is None:
		# return only cols equal below sig_cutoff
		sig_cutoff = 0.05

	if n_resamples is None:
		# return only cols equal below sig_cutoff
		n_resamples = 1000

	# Determine what/how many resamples need to be added to db
	toAdd = [ check_colResamples( i, len( rows ), n_resamples, host, port , db ) for i in cols]
	toAdd = [ i for i in toAdd if i is not None]

	count = 1
	if len( toAdd) > 0:
		if add_override:
			print "I need to perform %i random resample(s) of size %i to compute pvals. Please be patient. This may take a while..." % ( len(toAdd), n_resamples )
			tmp = colResampleInd( host = host, n_rows = len(rows), cols = toAdd, n_resamples = n_resamples, keepP = keepP, port = port, db = db )
			print "Done adding random resamples."
		else:
			print "I would need to perform %i random resample(s) of size %i to compute pvals. Since this would require significant computational power (and time), I have only returned results where resample data has been pre-calculated. Consult resample.py to run these jobs on multiple cores (much faster) or change 'add_override' flag of this function to 'True' to build the resample now." % ( len(toAdd), n_resamples )
			cols = [i for i in cols if i not in toAdd]

	if verbose:
		print "Calculating pvals"

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	exp_df = pd.DataFrame( list( client[ db ].gene_expression.find( { "col_id": { "$in" : cols }, "row_id": { "$in" : rows } }, { "_id":0, "col_id":1, "normalized_expression":1, "standardized_expression":1 } ) ) )
	random_rsd = pd.DataFrame( list( client[ db ].col_resample.find( { "n_rows": len( rows ), "col_id": { "$in" : cols } }, { "_id":0 } ) ) )
	
	if random_rsd.shape[0] == 0:
		print "Could not find resample DB entry for %i rows in cols %s" % ( len(rows_o), cols_o )
		return None
	else:
		random_rsd.index = random_rsd["col_id"]

	exp_df_rsd = exp_df.groupby("col_id").aggregate(rsd)

	if standardized:
		exp_df_rsd = exp_df_rsd.loc[ :, "standardized_expression"]
		resamples = random_rsd.loc[ :, "resamples"].to_dict()
		random_rsd = random_rsd.loc[ :,"lowest_standardized" ].to_dict()
	else:
		exp_df_rsd = exp_df_rsd.loc[ :, "raw_expression"]
		resamples = random_rsd.loc[ :, "resamples"].to_dict()
		random_rsd = random_rsd.loc[ :,"lowest_raw" ].to_dict()
		

	pvals = exp_df_rsd.groupby( level=0 ).aggregate(empirical_pval, random_rsd, resamples )
	pvals.columns = ["pval"]
	pvals.index = col2id_batch( pvals.index.values, host, port, db, input_type = "col_id", return_field = col_outtype )

	if sig_cutoff is not None:
		pvals = pvals[ pvals <= sig_cutoff ]

	if sort:
		pvals.sort()

	pvals = pvals.to_frame()
	pvals.columns = [ "pval" ] 

	if pvals.shape[ 0 ] == 0:
		print "No cols pass the significance cutoff of %f" % sig_cutoff

	return pvals

def agglom( x = [ 0,1 ], x_type = None, y_type = None, x_input_type = None, y_output_type = None, logic = "and", host = "localhost", port = 27017, db = "",  verbose = False, gre_lim = 10, pval_cutoff = 0.05, translate = True ):
	"""
	Determine enrichment of y given x through bicluster co-membership. 

	Available x_type(s)/y_type(s):

	'rows' or 'genes': search for row (genes) in biclusters. x should be a list of rows, eg ["carA","carB"] or [275, 276]
	'columns' or 'conditions': search for columns (conditions) in biclusters. x should be a list of columns, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
	'gre': search for GREs in biclusters. x should be a list of GRE IDs, eg [4, 19]
	'bicluster': self explanatory. takes or outputs bicluster '_id'
	''

	"""

	def compute_p( i, M, N ):
		#print i
		z = i.counts # n black balls in draw
		n = i.all_counts # num black balls tot
		M = M
		#M = M - n  # n white balls
		N = N # num drawn
		prb =  hypergeom.sf(z, M, n, N)
		return prb

	if x_type is None:
		print "Please supply an x_type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres' "
	if y_type is None:
		print "Please supply an y_type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres'. Biclusters will be returned by default. "
		y_type = "cluster"

	if type( x ) == str or type( x ) == int:
		# single
		x = [ x ]

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	# Check input types

	if x_type == "rows" or x_type == "row" or x_type == "gene" or x_type == "genes":
		x_type = "rows"
		x_o = x
		x = row2id_batch( x, host, port, db, input_type = x_input_type, return_field="row_id" )
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate row names: %s" % x_o
			return None 
	elif x_type == "columns" or x_type == "column" or x_type == "col" or x_type == "cols" or x_type == "condition" or x_type == "conditions" or x_type - "conds":
		x_type = "columns"
		x_o = x
		x = row2id_batch( x, host, port, db, input_type = x_input_type, return_field="col_id" )
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate col names: %s" % x_o
			return None 
	elif x_type == "motif" or x_type == "gre" or x_type == "motc" or x_type == "motif.gre" or x_type == "motifs" or x_type == "gres" or x_type == "motcs":
		x_type = "gre_id"
	elif x_type == "cluster" or x_type == "clusters" or x_type == "bicluster" or x_type == "biclusters" or x_type == "bcs":
		print "WARNING! I hope you are using cluster '_id'!!! Otherwise the results might surprise you..."
		x_type = "_id"
	else:
		print "ERROR: Can't recognize your 'x_type' argument."
		return None

	# Check output types

	if y_type == "rows" or y_type == "row" or y_type == "gene" or y_type == "genes":
		y_type = "rows"
	elif y_type == "columns" or y_type == "column" or y_type == "col" or y_type == "cols" or y_type == "condition" or y_type == "conditions" or y_type == "conds":
		y_type = "columns"
	elif y_type == "motif" or y_type == "gre" or y_type == "motc" or y_type == "motif.gre" or y_type == "motfs" or y_type == "gres" or y_type == "motcs":
		y_type = "gre_id"
	elif y_type == "cluster" or y_type == "clusters" or y_type == "bicluster" or y_type == "biclusters" or x_type == "bcs":
		print "WARNING! Will return bicluster _id. The results might surprise you..."
		y_type = "_id"
	else:
		print "ERROR: Can't recognize your 'y_type' argument."
		return None

	# Compose query

	if logic in [ "and","or","nor" ]:
		q = { "$"+logic: [ { x_type : i } for i in x ] }
		o = { y_type: 1 }
		if x_type == "gre_id":
			queryPre = pd.DataFrame( list( client[db].motif_info.find( q, { "cluster_id" : 1 } ) ) )["cluster_id"].tolist()
			query = pd.DataFrame( list( client[db].bicluster_info.find( { "_id": { "$in": queryPre } }, o ) ) )
		elif y_type == "gre_id":
			queryPre = pd.DataFrame( list( client[db].bicluster_info.find( q, { "_id" : 1 } ) ) )["_id"].tolist()
			query = pd.DataFrame( list( client[db].motif_info.find(  { "cluster_id": { "$in": queryPre } }, { y_type : 1 } ) ) )
		else:
			query = pd.DataFrame( list( client[db].bicluster_info.find( q, o ) ) )
	else:
		print "I don't recognize the logic you are trying to use. 'logic' must be 'and', 'or', or 'nor'."
		return None
	
	client.close()

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

				if client[db].rowsCount_mapreduce.count() == 0:
					print "Initializing MapReduce lookup table. Future queries will be much faster!"
					client[db].bicluster_info.map_reduce(mapRows,reduce,"rowsCount_mapreduce")
				else:
					# do spot check to make sure mapreduce is up to date
					random_id = random.randint( 0, client[db].rowsCount_mapreduce.count() )
					ref = client[db].rowsCount_mapreduce.find_one( { "_id": random_id } )["value"]
					test = client[db].bicluster_info.find( { "rows": random_id } ).count()
					if ref != test:
						print "Initializing MapReduce lookup table. Future queries will be much faster!"
						client[db].bicluster_info.map_reduce(mapRows,reduce,"rowsCount_mapreduce")

				rows = pd.Series( list( itertools.chain( *query.rows.tolist() ) ) ).value_counts().to_frame( "counts" )
				
				# filter out rows that aren't in the database - i.e. not annotated in MicrobesOnline
				in_db = pd.DataFrame( list( client[db].row_info.find( { }, { "_id" : 0, "row_id": 1 } ) ) ).row_id.tolist()
				common_rows = list(set(rows.index).intersection(set(in_db)))
				rows = rows.loc[common_rows]
				
				# find all bicluster counts
				all_counts = pd.DataFrame( list( client[db].rowsCount_mapreduce.find( ) ) )
				all_counts = all_counts.set_index("_id")

				# combine two data frames
				to_r = rows.join(all_counts).sort("counts",ascending=False)
				to_r.columns = ["counts","all_counts"]

				if translate:
					to_r.index = row2id_batch( to_r.index.tolist(), host, port, db, return_field = "egrin2_row_name", input_type = "row_id" )

			if y_type == "columns":

				if client[db].columnsCount_mapreduce.count() == 0:
					print "Initializing MapReduce lookup table. Future queries will be much faster!"
					client[db].bicluster_info.map_reduce(mapColumns,reduce,"columnsCount_mapreduce")
				else:
					# do spot check to make sure mapreduce is up to date
					random_id = random.randint( 0, client[db].columnsCount_mapreduce.count() )
					ref = client[db].columnsCount_mapreduce.find_one( { "_id": random_id } )["value"]
					test = client[db].bicluster_info.find( { "columns": random_id } ).count()
					if ref != test:
						print "Initializing MapReduce lookup table. Future queries will be much faster!"
						client[db].bicluster_info.map_reduce(mapColumns,reduce,"columnsCount_mapreduce")

				if client[db].columnsCount_mapreduce.count() == 0:
					client[db].bicluster_info.map_reduce(mapColumns,reduce,"columnsCount_mapreduce")
				
				cols = pd.Series( list( itertools.chain( *query["columns"].tolist() ) ) ).value_counts().to_frame( "counts" )
				
				# filter out columns that aren't in the database - i.e. not annotated in MicrobesOnline
				in_db = pd.DataFrame( list( client[db].col_info.find( { }, { "_id" : 0, "col_id": 1 } ) ) ).col_id.tolist()
				common_cols = list( set( cols.index ).intersection( set( in_db )  ))
				cols = cols.loc[common_cols]
				
				# find all bicluster counts
				all_counts = pd.DataFrame( list( client[db].columnsCount_mapreduce.find( ) ) )
				all_counts = all_counts.set_index("_id")

				# combine two data frames
				to_r = cols.join(all_counts).sort("counts",ascending=False)
				to_r.columns = ["counts","all_counts"]

				if translate:
					to_r.index = col2id_batch( to_r.index.tolist(), host, port, db, return_field = "egrin2_col_name", input_type = "col_id" )


			if y_type == "gre_id":

				if client[db].gresCount_mapreduce.count() == 0:
					print "Initializing MapReduce lookup table. Future queries will be much faster!"
					client[db].motif_info.map_reduce(mapGREs,reduce,"gresCount_mapreduce")
				else:
					# do spot check to make sure mapreduce is up to date
					random_id = random.randint( 0, client[db].gresCount_mapreduce.count() )
					ref = client[db].gresCount_mapreduce.find_one( { "_id": random_id } )["value"]
					test = client[db].motif_info.find( { "gre_id": random_id } ).count()
					if ref != test:
						print "Initializing MapReduce lookup table. Future queries will be much faster!"
						client[db].motif_info.map_reduce(mapGREs,reduce,"gresCount_mapreduce")

				gres = query.gre_id.tolist() 
				gres = filter(lambda x: x != "NaN", gres )
				gres = pd.Series( gres ).value_counts().to_frame( "counts" )

				# find all bicluster counts
				all_counts = pd.DataFrame( list( client[db].gresCount_mapreduce.find( ) ) )
				all_counts = all_counts.set_index("_id")

				# combine two data frames
				to_r = gres.join(all_counts).sort("counts",ascending=False)
				to_r.columns = ["counts","all_counts"]

				# filter by GREs with more than 10 instances
				to_r = to_r.loc[ to_r.all_counts>=gre_lim, : ]


			to_r["pval"] = to_r.apply( compute_p, axis=1, M = client[db].bicluster_info.count(),  N = query.shape[ 0 ] )
			to_r["qval_BH"] = multipletests( to_r.pval, method='fdr_bh' )[1]
			to_r["qval_bonferroni"] = multipletests( to_r.pval, method='bonferroni' )[1]
			to_r = to_r.sort( "pval", ascending=True )
			# only return below pval cutoff
			to_r = to_r.loc[ to_r.pval <= pval_cutoff, : ]

			return to_r
	
	else:
		
		print "Could not find any biclusters matching your criteria"
		return None

def fimoFinder( start = None, stop = None, chr = None, strand = None, mot_pval_cutoff = None, filter = None, filter_type = None, host = "localhost", port = 27017, db = "" ):
	"""Find motifs/GREs that fall within a specific range. Filter by biclusters/genes/conditions/etc."""
	
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	mots = pd.DataFrame( list( client[db].fimo.find( { "start": { "gte": start }, "stop": {"lte": stop } } ) ) )

def coremFinder( x, x_type = "corem_id", x_input_type = None, y_type = "genes", y_return_field = None, count = False, logic = "or", host = "localhost", port = 27017, db = "" ):

	if x_type is None:
		print "Please supply an x_type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres', edges "
	if y_type is None:
		print "Please supply an y_type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres', edges "

	if type( x ) == str or type( x ) == int:
		# single
		x = [ x ]

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	# Check input types

	if x_type == "rows" or x_type == "row" or x_type == "gene" or x_type == "genes":
		x_type = "rows"
		x_o = x
		x = row2id_batch( x, host, port, db, input_type = x_input_type, return_field="row_id" )
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate row names: %s" % x_o
			return None 
	elif x_type == "columns" or x_type == "column" or x_type == "col" or x_type == "cols" or x_type == "condition" or x_type == "conditions" or x_type == "conds":
		x_type = "cols"
		x_o = x
		x = col2id_batch( x, host, port, db, input_type = x_input_type, return_field="col_id" )
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate row names: %s" % x_o
			return None 
	elif x_type == "motif" or x_type == "gre" or x_type == "motc" or x_type == "motif.gre" or x_type == "motifs" or x_type == "gres" or x_type == "motcs":
		x_type = "gre_id"
	elif x_type == "corem_id" or x_type == "corem" or "corems":
		x_type = "corem_id"
	
	# Check output types

	if y_type == "rows" or y_type == "row" or y_type == "gene" or y_type == "genes":
		y_type = "rows"

	elif y_type == "columns" or y_type == "column" or y_type == "col" or y_type == "cols" or y_type == "condition" or y_type == "conditions" or y_type == "conds":
		y_type = "cols"

	if logic in [ "and","or","nor" ]:
		if logic == "and" and x_type == "corem_id":
			q = { "$or": [ { x_type : i } for i in x ] }
		else:
			q = { "$"+logic: [ { x_type : i } for i in x ] }
		o = { y_type: 1 }
		query = pd.DataFrame( list( client[db].corem.find( q, o ) ) )
	else:
		print "I don't recognize the logic you are trying to use. 'logic' must be 'and', 'or', or 'nor'."
		return None

	if y_type == "rows":
		if y_return_field is None:
			y_return_field = "egrin2_row_name"
		if query.shape[0] > 1:
			to_r = list( itertools.chain( *query.rows.values.tolist() ) )
			if logic == "and":
				to_r = pd.Series( to_r ).value_counts()
				if count:
					to_r = to_r[ to_r >= query.shape[0] ]
					if to_r.shape[0] > 0:
						to_r.index = row2id_batch( to_r.index.tolist(), host, port, db, return_field = y_return_field, input_type = "row_id" )
					else:
						print "No genes found"
						return None
				else:
					to_r = to_r[ to_r > query.shape[0] ].index.tolist()
					if len( to_r ):
						to_r = row2id_batch( to_r, host, port, db, return_field = y_return_field, input_type = "row_id" )
						to_r.sort()
					else:
						print "No genes found"
						return None
			else:
				if count:
					to_r = pd.Series( to_r ).value_counts()
					if to_r.shape[0] > 0:
						to_r.index = row2id_batch( to_r.index.tolist(), host, port, db, return_field = y_return_field, input_type = "row_id" )
					else:
						print "No genes found"
						return None
				else:
					to_r = list( set( to_r ) )
					if len( to_r ):
						to_r = row2id_batch( to_r, host, port, db, return_field = y_return_field, input_type = "row_id" )
						to_r.sort()
					else:
						print "No genes found"
						return None 
		else:
			to_r = row2id_batch( query.rows[0], host, port, db, return_field = y_return_field, input_type = "row_id" )

	elif y_type == "cols":
		if y_return_field is None:
			y_return_field = "egrin2_col_name"
		if query.shape[0] > 1:
			to_r = [int(i["col_id"]) for i in list(itertools.chain( *query.cols.values.tolist())) if i["col_id"] if type(i["col_id"]) is float]
			if logic == "and":
				to_r = pd.Series( to_r ).value_counts()
				if count:
					to_r = to_r[ to_r >= query.shape[0] ]
					if to_r.shape[0] > 0:
						to_r.index = col2id_batch( to_r.index.tolist(), host, port, db, return_field = y_return_field, input_type = "col_id" )
					else:
						print "No conditions found"
						return None
				else:
					to_r = to_r[ to_r >= query.shape[0] ].index.tolist()
					if len( to_r ):
						to_r = col2id_batch( to_r, host, port, db, return_field = y_return_field, input_type = "col_id" )
						to_r.sort()
					else:
						print "No conditions found"
						return None
			else:
				if count:
					to_r = pd.Series( to_r ).value_counts()
					if to_r.shape[0] > 0:
						to_r.index = col2id_batch( to_r.index.tolist(), host, port, db, return_field = y_return_field, input_type = "col_id" )
					else:
						print "No conditions found"
						return None
				else:
					to_r = list( set( to_r ) )
					if len( to_r ):
						to_r = col2id_batch( to_r, host, port, db, return_field = y_return_field, input_type = "col_id" )
						to_r.sort()
					else:
						print "No conditions found"
						return None 
		else:
			to_r = [int(i["col_id"]) for i in list(itertools.chain( *query.cols.values.tolist())) if i["col_id"] if type(i["col_id"]) is float]
			to_r = col2id_batch( to_r, host, port, db, return_field = y_return_field, input_type = "col_id" )
	else:
		print "Could not find corems matching your query"
		to_r = None

	return to_r

def expressionFinder( rows = None, cols = None, standardized = True, host = "localhost", port = 27017, db = "" ):
	"""Find motifs/GREs that fall within a specific range. Filter by biclusters/genes/conditions/etc."""
	
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	input_type_rows = None
	input_type_cols = None
	if rows is None:
		# assume all genes
		rows = pd.DataFrame( list( client[ db ].row_info.find( {}, {"row_id":1} ) ) ).row_id.tolist()
		input_type_rows = "row_id"
	if cols is None:
		# assume all cols
		cols = pd.DataFrame( list( client[ db ].col_info.find( {}, {"col_id":1} ) ) ).col_id.tolist()
		input_type_cols = "col_id"

	if type( rows ) == str or type( rows ) == int:
		# single
		rows = [ rows ]

	if type( cols ) == str or type( cols ) == int:
		# single
		cols = [ cols ]

	# translate rows/cols

	rows = row2id_batch( rows, host, port, db,  verbose = False, return_field = "row_id", input_type = input_type_rows )
	cols = col2id_batch( cols, host, port, db,  verbose = False, return_field = "col_id", input_type = input_type_cols )

	if len( rows ) > 1000 or len( cols ) > 1000:
		print "WARNING: This is a large query. Please be patient. If you need faster access, I would suggest saving this matrix and loading directly from file."

	# get expression data
	data = pd.DataFrame( None,columns = cols, index = rows )
	query = pd.DataFrame ( list( client[ db ].gene_expression.find( { "$and": [ { "row_id": { "$in": rows } }, { "col_id": { "$in": cols } } ] } ) ) )
	for i in query:
		if standardized:
			data = query.pivot(index="row_id",columns="col_id",values="standardized_expression")
		else:
			data = query.pivot(index="row_id",columns="col_id",values="raw_expression")

	data.index = row2id_batch( data.index.tolist(), host, port, db,  verbose = False, return_field = "egrin2_row_name", input_type = "row_id" )
	data.columns = col2id_batch( data.columns.tolist(), host, port, db,  verbose = False, return_field = "egrin2_col_name", input_type = "col_id" )
	data = data.sort_index( )
	data = data.reindex_axis( sorted( data.columns ), axis=1 )

	return data

def plotExpression( data, plot_type = "boxplot", ipynb = False, zlim = None, sort = True, boxpoints = None ):
	
	if sort:
		c_order = data.mean(0).order().index.tolist()
		data = data.loc[:,c_order]

	def to_scatter( df ):
		x = df.index.values
		lines={}
		for key in df:
		    lines[key]={}
		    lines[key]["x"]=x
		    lines[key]["y"]=df[key].values
		    lines[key]["name"]=key

		    #Appending all lines
		lines_plotly=[lines[key] for key in df]
		return lines_plotly

	def to_box( df, boxpoints ):

		if boxpoints:
			boxpoints = "all"
		elif df.shape[1]<=50 and boxpoints is None:
			boxpoints = "all"
		elif boxpoints is None:
			boxpoints = False
		else:
			boxpoints = boxpoints

		boxes = []
		for x in df.columns.tolist():
			boxes.append(Box(
					name = x,
					y = df.loc[:,x].tolist(),
					boxpoints=boxpoints,
        			jitter=0.3,
        			pointpos=0
				))
		return boxes

	def to_heatmap( df, zlim ):

		if zlim is None:
			# use 1st and 99th percentile
			all_data = list( itertools.chain( *df.values.tolist() ) )
			if np.percentile( all_data, 1 ) < 0:
				zmin = -max( np.abs( np.percentile( all_data, [ 1, 99 ] ) ) )
			else:
				zmin = np.percentile( all_data, 1 )
			zmax = max( np.abs( np.percentile( all_data, [ 1, 99 ] ) ) )
		else:
			zmin = zlim[0]
			zmax = zlim[1]

		# do hierarchical clustering to make it look pretty
		D1 = squareform(pdist(df, metric='euclidean'))
		D2 = squareform(pdist(df.T, metric='euclidean'))
		Y = linkage(D1, method='complete')
		Z1 = dendrogram(Y, orientation='right')
		Y = linkage(D2, method='complete')
		Z2 = dendrogram(Y)
		idx1 = Z1['leaves']
		idx2 = Z2['leaves']
		D = df.iloc[idx1, :]
		D = D.iloc[:, idx2]
		to_r = Heatmap(
		        z = [ D.loc[i].tolist() for i in D.index.tolist() ],
		        y = D.index.tolist(),
		        x = D.columns.tolist(),
		        colorscale=[ [0.0, 'rgb(0,0,255)'], [0.5, 'rgb(0,0,0)'], [1.0, 'rgb(255,255,0)'] ],
		        zauto = False,
		        zmin = zmin,
		        zmax = zmax

		    )
		return to_r

	if plot_type == "line":
		to_plot = to_scatter(data.T)

		layout = Layout(
		title='Expression',
		xaxis=XAxis(
			title='Condition',
			ticks='',
			showticklabels=False
			),
		yaxis=YAxis(
			title='Expression Value',
			zeroline=True,
			)
		)

		fig = Figure(data=to_plot, layout=layout)

	elif plot_type == "heatmap":
		to_plot = to_heatmap(data, zlim)
		fig = Data( [ to_plot ] )
	elif plot_type == "boxplot":
		to_plot = to_box(data, boxpoints)
		fig = to_plot
	else:
		print "ERROR: Cannot recognize plot_type = %s" % plot_type
		return None


	if not ipynb:
		unique_url = py.plot( fig )
		return unique_url

	return fig






if __name__ == '__main__':
	print "yuuuup"
	# host = "primordial"
	# port = 27017
	# db = "egrin2_db"

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	# test rows
	corem = client[db].corem.find_one({"corem_id":1})["rows"]
	# bicluster_ids = x2bicluster( x = corem, x_type = "rows", logic = "or", count = False, host = host, port = port, db = db,  verbose = False, return_field = [ "_id" ] )._id.tolist()

	# # test cols
	# cols = range(5)
	# bicluster_ids = x2bicluster( x = cols, x_type = "cols", logic = "or", count = False, host = host, port = port, db = db,  verbose = False, return_field = [ "_id" ] )._id.tolist()

	# # test gres
	# gre = [2]
	# bicluster_ids = x2bicluster( x = gre, x_type = "gre", logic = "or", count = False, host = host, port = port, db = db,  verbose = False, return_field = [ "_id" ] )._id.tolist()