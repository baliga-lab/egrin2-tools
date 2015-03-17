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

from pymongo import MongoClient
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import itertools
from bson.code import Code
import matplotlib.pyplot as plt
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
		if return_field == "all":
			row = query[0]
		else:
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
		to_r= [ row2id( x, host, port, db, return_field = return_field ) for x in rows ]

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
		to_r= [ col2id( x, host, port, db, return_field = return_field ) for x in cols ]
	
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

	client.close()

	return pvals

def agglom( x = [ 0,1 ], x_type = None, y_type = None, x_input_type = None, y_output_type = None, logic = "or", host = "localhost", port = 27017, db = "",  verbose = False, gre_lim = 10, pval_cutoff = 0.05, translate = True ):
	"""
	Determine enrichment of y given x through bicluster co-membership. 

	Available x_type(s)/y_type(s):

	'rows' or 'genes': x should be a gene name or list of genes, eg ["carA","carB"] or [275, 276]
	
	'columns' or 'conditions': x should be a condition name or list of conditions, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
	
	'gre':  x should be a GRE ID or list of GRE IDs, eg [4, 19]
	
	'bicluster': takes or outputs bicluster '_id'
	''

	"""

	print 'Using "%s" logic' % logic

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
	elif x_type == "columns" or x_type == "column" or x_type == "col" or x_type == "cols" or x_type == "condition" or x_type == "conditions" or x_type == "conds":
		x_type = "columns"
		x_o = x
		x = col2id_batch( x, host, port, db, input_type = x_input_type, return_field="col_id" )
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
			to_r = to_r.sort( ["pval","counts"], ascending=True )
			# only return below pval cutoff
			to_r = to_r.loc[ to_r.pval <= pval_cutoff, : ]

			return to_r
	
	else:
		
		print "Could not find any biclusters matching your criteria"
		return None

def fimoFinder( start = None, stop = None, locusId = None, strand = None, mot_pval_cutoff = None, filterby = None, filter_type = None, filterby_input_type = None, host = "localhost", port = 27017, db = None, use_fimo_small = True, logic = "or", return_format = "file", outfile = None, tosingle = True ):
	"""Find motifs/GREs that fall within a specific range. Filter by biclusters/genes/conditions/etc."""
	
	def getBCs( x, x_type ):
		if x is None:
			to_r = pd.DataFrame( list( client[ db ].motif_info.find( { }, { "cluster_id" : 1, "gre_id" : 1} ) ) )
		else:
			to_r = pd.DataFrame( list( client[ db ].motif_info.find( { x_type: x }, { "cluster_id" : 1, "gre_id" : 1 } ) ) )
		return( to_r.loc[ :, [ "gre_id" ,"cluster_id" ] ] )

	def aggSeq( x ):
		# print x.gre_id
		def count( y ):
			return( range( y.start, y.stop+1) )
		if x.shape[ 0 ] == 0:
			return( [] )
		else:
			to_r = [ count( x.iloc[ i ] ) for i in range( x.shape[ 0 ] ) ]
			to_r = list( itertools.chain( *to_r ) )
			to_r.sort()
		return( to_r )

	try:
		client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	except Exception:
		print "Cant connect to MongoDB at host = %s, port = %s" % ( host, str(port) )

	if db is None:
		print "Please provide a database name, e.g. *org*_db, where *org* is a three lettter short organism code"
		return None

	db_chr = pd.DataFrame( list( client[ db ].genome.find( { }, { "scaffoldId":1, "NCBI_RefSeq":1 } ) ) )
	db_scaffoldId = db_chr.scaffoldId.tolist( )
	db_NCBI_RefSeq = db_chr.NCBI_RefSeq.tolist( )
	
	if locusId is None:
		print "Please provide a chromosome Locus ID. This is probably the scaffoldID from MicrobesOnline.\n\nLocusIds in database %s include: \n\nScaffoldId\n%s\n\nNCBI_RefSeq\n%s" % ( db, (", ").join( db_scaffoldId ), (", ").join( db_NCBI_RefSeq ) )
		return None

	locusId = str( locusId )
	
	if locusId not in db_scaffoldId and locusId not in db_NCBI_RefSeq:
		print "LocusId %s not in EGRIN 2.0 database %s. \n\nLocusIds in this database include: \n\nScaffoldId\n%s\n\nNCBI_RefSeq\n%s" % ( locusId, db, (", ").join( db_scaffoldId ), (", ").join( db_NCBI_RefSeq ) )
		return None

	chromosome = client[ db ].genome.find_one( { "$or": [ { "scaffoldId": locusId }, { "NCBI_RefSeq": locusId } ] } )
	scaffoldId = chromosome[ "scaffoldId" ]
	ncbi = chromosome[ "NCBI_RefSeq" ]

	# if start/stop is None, assume whole chromosome
	if start is None:
		print "Start not provided. Assuming beginning of chromosome"
		start = 0

	if stop is None:
		print "Stop not provided. Assuming end of chromosome"
		stop = len( chromosome["sequence"] )

	if use_fimo_small:
		fimo_collection = "fimo_small"
	else:
		fimo_collection = "fimo"

	if filterby is None:
		print "No filter applied"
	
	if filter_type is not None:
		print "WARNING: Many of these filters are not supported currently. Only GREs!!!"
		if filter_type == "rows" or filter_type == "row" or filter_type == "gene" or filter_type == "genes":
			filter_type = "rows"
			filterby_o = filterby
			filterby = row2id_batch( filterby, host, port, db, input_type = filter_input_type, return_field="row_id" )
			filterby = list( set( filterby ) )
			if len( filterby ) == 0:
				print "Cannot translate row names: %s" % filterby_o
				return None 
		elif filter_type == "columns" or filter_type == "column" or filter_type == "col" or filter_type == "cols" or filter_type == "condition" or filter_type == "conditions" or filter_type == "conds":
			filter_type = "columns"
			filterby_o = filterby
			filterby = col2id_batch( filterby, host, port, db, input_type = filterby_input_type, return_field="col_id" )
			filterby = list( set( filterby ) )
			if len( filterby ) == 0:
				print "Cannot translate col names: %s" % filterby_o
				return None 
		elif filter_type == "motif" or filter_type == "gre" or filter_type == "motc" or filter_type == "motif.gre" or filter_type == "motifs" or filter_type == "gres" or filter_type == "motcs":
			filter_type = "gre_id"
		elif filter_type == "cluster" or filter_type == "clusters" or filter_type == "bicluster" or filter_type == "biclusters" or filter_type == "bcs":
			print "WARNING! I hope you are using cluster '_id'!!! Otherwise the results might surprise you..."
			filter_type = "_id"
		print "Filtering motifs by %s" % filter_type
		bcs_df= pd.concat( [ getBCs( i, filter_type ) for i in filterby ], ignore_index = True ) 

		mots = pd.DataFrame( list( client[db][fimo_collection].find( { "start": { "$gte": start }, "stop": {"$lte": stop }, "cluster_id": {"$in": bcs_df.cluster_id.tolist() }, "scaffoldId": scaffoldId } ) ) )
		mots = pd.merge( mots, bcs_df, on= "cluster_id" )
	else:
		bcs_df= pd.concat( [ getBCs( None, filter_type ) ], ignore_index = True ) 
		mots = pd.DataFrame( list( client[db][fimo_collection].find( { "start": { "$gte": start }, "stop": {"$lte": stop }, "scaffoldId": scaffoldId } ) ) )
		mots = pd.merge( mots, bcs_df, on= "cluster_id" )

	gre_scans = mots.groupby( "gre_id" ).apply( aggSeq )

	if return_format == "dictionary":
		gre_scans = { i: pd.Series( gre_scans.loc[ i ] ).value_counts().sort_index() for i in gre_scans.index }

	if return_format == "file":
		# return with file format ready to save for GGBweb
		# start stop strand chr value id
		tmp_gre_scans = []
		for i in gre_scans.index:
			tmp_df = pd.DataFrame(columns=[ "start","end","strand","chr","value","id" ] )
			gre_counts = pd.Series( gre_scans.loc[ i ] ).value_counts().sort_index()
			tmp_df[ "start" ] =  gre_counts.index
			tmp_df[ "end" ] =  gre_counts.index
			tmp_df[ "strand" ] = "+"
			tmp_df[ "chr" ] = ncbi
			tmp_df[ "value" ] = gre_counts.values
			tmp_df[ "id" ] = "GRE_" + str(i)
			tmp_gre_scans.append( tmp_df )
		gre_scans = pd.concat( tmp_gre_scans, ignore_index=True )

		if outfile is not None:
			def dfsave( df, fname ):
				fname = list( set( df.id ) )[ 0 ]  + "_" + fname
				print "Writing file %s" % fname
				df = df.drop('id', 1)
				df.to_csv( fname, sep="\t", index=False )
				return None
			if tosingle:
				gre_scans_grouped = gre_scans.groupby( "id" )
				tmp = [ dfsave( gre_scans_grouped.get_group( i ), fname=outfile ) for i in gre_scans_grouped.groups.keys( ) ]
				return None
			else:
				gres_scans.to_csv(  outfile, sep="\t", index=False )
				return None

	client.close()

	return gre_scans

def coremFinder( x, x_type = "corem_id", x_input_type = None, y_type = "genes", y_return_field = None, count = False, logic = "or", host = "localhost", port = 27017, db = "" ):

	"""
	Fetch corem-related info 'y' given query 'x'. 

	Available x_type(s)/y_type(s):

	'corem_id': x is corem_id (int), eg [1]

	'rows' or 'genes': x should be a gene name or list of genes, eg ["carA","carB"] or [275, 276]
	
	'columns' or 'conditions': x should be a condition or list of conditions, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
	
	'gre': x should be a GRE or a list of GRE IDs, eg [4, 19]

	'edge': x should be an edge or list of edges. Genes in edges should be separated by '-', eg ["carA-carB"] or ["275-276"]
	
	"""


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
		x_type = "cols.col_id"
		x_o = x
		x = col2id_batch( x, host, port, db, input_type = x_input_type, return_field="col_id" )
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate row names: %s" % x_o
			return None 
	elif x_type == "motif" or x_type == "gre" or x_type == "motc" or x_type == "motif.gre" or x_type == "motifs" or x_type == "gres" or x_type == "motcs":
		x_type = "gre_id"
	elif x_type == "corem_id" or x_type == "corem" or x_type == "corems":
		x_type = "corem_id"
	elif x_type == "edge" or x_type == "edges":
		x_type = "edges"
		x_new = []
		for i in x:
			i_trans = row2id_batch( i.split("-"), host, port, db, input_type = x_input_type, return_field="row_id", verbose = False )
			i_trans = [ str( j ) for j in i_trans ]
			x_new.append( "-".join( i_trans ) )
			i_trans.reverse()
			x_new.append( "-".join( i_trans ) )
		x = x_new
		if len( x ) == 0:
			print "Cannot translate row names: %s" % x_o
			return None 

	# Check output types

	if y_type == "rows" or y_type == "row" or y_type == "gene" or y_type == "genes":
		y_type = "rows"
	elif y_type == "columns" or y_type == "column" or y_type == "col" or y_type == "cols" or y_type == "condition" or y_type == "conditions" or y_type == "conds":
		y_type = "cols.col_id"
	elif y_type == "motif" or y_type == "gre" or y_type == "motc" or y_type == "motif.gre" or y_type == "motifs" or y_type == "gres" or y_type == "motcs":
		y_type = "gre_id"
	elif y_type == "corem_id" or y_type == "corem" or y_type == "corems":
		y_type = "corem_id"
	elif y_type == "edge" or y_type == "edges":
		y_type = "edges"

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

	elif y_type == "cols.col_id":
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
	elif y_type == "corem_id":
		to_r = query.corem_id.tolist()
	elif y_type == "edges":
		if y_return_field is None:
			y_return_field = "egrin2_row_name"
		to_r = list( itertools.chain( *query.edges.values.tolist() ) )
		to_r_new = []
		for i in to_r:
			i_trans = row2id_batch( [ int( j ) for j in i.split("-") ], host, port, db, input_type = "row_id", return_field=y_return_field, verbose = False )
			i_trans = [ str( j ) for j in i_trans ]
			to_r_new.append( "-".join( i_trans ) )
		to_r = to_r_new
		to_r.sort()
	elif y_type == "gre_id":
		print "GREs detection is not currently supported for corems. Please use the `agglom` function to find GREs enriched in biclusters containing corem genes instead."
		to_r = None
	else:
		print "Could not find corems matching your query"
		to_r = None

	client.close()

	if to_r is not None:
		to_r = pd.DataFrame(to_r)

	return to_r

def expressionFinder( rows = None, cols = None, standardized = True, host = "localhost", port = 27017, db = "" ):
	"""Fetch gene expression given rows and columns."""
	
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

	client.close()

	return data

def ggbwebModule( genes = None, outfile = None, host = "localhost", port = 27017, db = "" ):
	
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	if genes is None:
		print "Please provide a gene or list of genes"
		return None

	if type( genes ) == str or type( genes ) == int:
		# single
		genes = [ genes ]

	def locFormat(x):
    		return( "chromosome" + x.strand + ":" + str( x.start ) + "-" + str( x.stop ) )
	
	gene_info = pd.DataFrame( row2id_batch( genes, host=host, port=port, db=db,  return_field = "all", verbose = False ) )

	to_r = pd.DataFrame( gene_info.loc[ :, "egrin2_row_name" ] )
	to_r[ "loc" ] = [ locFormat( gene_info.iloc[ i ] ) for i in range( gene_info.shape[ 0 ] ) ]
	to_r[ "name2" ] = gene_info[ "name" ]

	if outfile is not None:
		print "Module written to: %s" % os.path.abspath( outfile )
		to_r.to_csv( os.path.abspath( outfile ), sep="\t", index=False, header=False)
		return None

	client.close()

	return to_r

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