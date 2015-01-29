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

from resample import *

def rsd( vals ):
	return abs( np.std( vals ) / np.mean( vals ) )

def check_colResamples( col, n_rows, n_resamples, host, port, db,  ):
	# connect to db
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' ) 
	if client[db].col_resample.find_one( { "n_rows": n_rows, "col_id": col, "resamples": { "$gte": n_resamples } } ) is None:
		client.close()
		return col
	client.close()

def row2id( row, host, port, db,  verbose = False, return_field = "row_id" ):
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

def row2id_batch( rows, host, port, db,  verbose = False, return_field = "row_id" ):
	"""Check name format of rows. If necessary, translate."""
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = pd.DataFrame( list( client[db].row_info.find( { "$or": [ { "row_id": { "$in": rows } }, { "egrin2_row_name": { "$in": rows } }, { "GI": { "$in": rows } }, { "accession": {"$in": rows } }, { "name": { "$in": rows } }, { "sysName": { "$in": rows }  } ] }, { return_field: 1 } ) ) ).loc[ :, return_field ].tolist()
	client.close()
	if len( rows ) > len( query ): 
		print "WARNING: Returning fewer rows than originally supplied"
	return query
	

def col2id( col, host, port, db,  verbose = False, return_field = "row_id" ):
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


def col2id_batch( cols, host, port, db,  verbose = False, return_field = "col_id" ):
	"""Check name format of rows. If necessary, translate."""
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = pd.DataFrame( list( client[db].col_info.find( { "$or": [ { "col_id": { "$in": cols } }, { "egrin2_col_name": { "$in":cols } } ] }, { return_field: 1 } ) ) ).loc[ :, return_field ].tolist()
	client.close()
	if len( cols ) > len( query ): 
		print "WARNING: Returning fewer rows than originally supplied"
	return query

def colResamplePval( rows = None, cols = None, n_resamples = None, host = "localhost", port = 27017, db = "", standardized = None, sig_cutoff = None, sort = True, add_override = False, n_jobs = 4, keepP = 0.1, verbose = True ):

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
	rows = [ row2id( i, host, port, db ) for i in rows ]
	rows = [ i for i in rows if i is not None]
	rows = list( set( rows ) )
	if len( rows ) == 0:
		print "Please provide an appropriately named array of rows"
		return None 
	
	cols_o = cols
	cols = [ col2id( i, host, port, db ) for i in cols ]
	cols = [ i for i in cols if i is not None]
	rows = list( set( rows ) )
	if len( rows ) == 0:
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
		exp_df_rsd = exp_df_rsd.loc[ :, "normalized_expression"]
		resamples = random_rsd.loc[ :, "resamples"].to_dict()
		random_rsd = random_rsd.loc[ :,"lowest_normalized" ].to_dict()
		

	pvals = exp_df_rsd.groupby( level=0 ).aggregate(empirical_pval, random_rsd, resamples )
	pvals.columns = ["pval"]
	pvals.index = [ col2name( i, host, port, db ) for i in pvals.index.values]

	if sig_cutoff is not None:
		pvals = pvals[ pvals <= sig_cutoff ]

	if sort:
		pvals.sort()

	pvals = pvals.to_frame()
	pvals.columns = [ "pval" ] 

	if pvals.shape[ 0 ] == 0:
		print "No cols pass the significance cutoff of %f" % sig_cutoff

	return pvals

def rows2corem( rows = [ 0, 1 ], host = "localhost", port = 27017, db = "",  verbose = False, return_field = [ "corem_id" ], logic = "and" ):
	"""Find corems in which row(s) co-occur."""
	
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	rows_o = rows
	rows = [ row2id( i, host, port, db ) for i in rows ]
	rows = [ i for i in rows if i is not None]
	rows = list( set( rows ) )
	if len( rows ) == 0:
		print "Cannot translate row names: %s" % rows_o
		return None 

	if logic in [ "and","or","nor" ]:
		q = { "$"+logic: [ { "rows" : i } for i in rows ] }
		query = pd.DataFrame( list( client[db].corem.find( q ) ) )
	else:
		print "I don't recognize the logic you are trying to use. 'logic' must be 'and', 'or', or 'nor'."
	
	client.close()

	if query.shape[0] > 0: 
		if return_field == "all":
			return query
		else:
			try:
				return query.loc[ :, return_field ]
			except Exception:
				return query
	else:
		print "Could not find any corems matching your criteria"
		return None

def agglom( x = [ 0,1 ], x_type = None, y_type = None, logic = "and", host = "localhost", port = 27017, db = "",  verbose = False, gre_lim = 10, pval_cutoff = 0.05, translate = True ):
	"""
	Determine how often 'x' occurs in biclusters. Usually just retrieve the counts. Retrieve additional bicluster info by setting count to False

	Available x_type(s)/y_type(s):

	'rows' or 'genes': search for row (genes) in biclusters. x should be a list of rows, eg ["carA","carB"] or [275, 276]
	'columns' or 'conditions': search for columns (conditions) in biclusters. x should be a list of columns, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
	'gre': search for GREs in biclusters. x should be a list of GRE IDs, eg [4, 19]
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
		x = [ row2id( i, host, port, db ) for i in x ]
		x = [ i for i in x if i is not None]
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate row names: %s" % x_o
			return None

	if x_type == "columns" or x_type == "column" or x_type == "col" or x_type == "cols" or x_type == "condition" or x_type == "conditions":
		x_type = "columns"
		x_o = x
		x = [ col2id( i, host, port, db ) for i in x ]
		x = [ i for i in x if i is not None]
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate col names: %s" % x_o
			return None  

	if x_type == "motif" or x_type == "gre" or x_type == "motc" or x_type == "motif.gre" or x_type == "motfs" or x_type == "gres" or x_type == "motcs":
		x_type = "motif.gre_id"

	if x_type == "cluster" or x_type == "clusters" or x_type == "bicluster" or x_type == "biclusters":
		print "WARNING! I hope you are using cluster _id!!! Otherwise the results might surprise you..."
		x_type = "_id"

	# Check output types

	if y_type == "rows" or y_type == "row" or y_type == "gene" or y_type == "genes":
		y_type = "rows"

	if y_type == "columns" or y_type == "column" or y_type == "col" or y_type == "cols" or y_type == "condition" or y_type == "conditions":
		y_type = "columns"
		
	if y_type == "motif" or y_type == "gre" or y_type == "motc" or y_type == "motif.gre" or y_type == "motfs" or y_type == "gres" or y_type == "motcs":
		y_type = "motif.gre_id"

	if y_type == "cluster" or y_type == "clusters" or y_type == "bicluster" or y_type == "biclusters":
		print "WARNING! Will return bicluster _id. The results might surprise you..."
		y_type = "_id"

	# Compose query

	if logic in [ "and","or","nor" ]:
		q = { "$"+logic: [ { x_type : i } for i in x ] }
		o = { y_type: 1 }
		query = pd.DataFrame( list( client[db].bicluster_info.find( q, o ) ) )
	else:
		print "I don't recognize the logic you are trying to use. 'logic' must be 'and', 'or', or 'nor'."
		return None
	
	client.close()

	if query.shape[0] > 0: 
		if y_type == "rows":
			rows = pd.Series( list( itertools.chain( *query.rows.tolist() ) ) ).value_counts().to_frame( "counts" )
			
			# filter out rows that aren't in the database - i.e. not annotated in MicrobesOnline
			in_db = pd.DataFrame( list( client[db].row_info.find( { }, { "_id" : 0, "row_id": 1 } ) ) ).row_id.tolist()
			common_rows = list(set(rows.index).intersection(set(in_db)))
			rows = rows.loc[common_rows]
			
			# find all bicluster counts
			all_counts = pd.DataFrame( list( client[db].bicluster_info.find( { }, { "rows": 1 } ) ) )
			all_counts = pd.Series( list( itertools.chain( *all_counts.rows.tolist() ) ) ).value_counts().to_frame( "all_counts" )
			# combine two data frames
			to_r = rows.join(all_counts).sort("counts",ascending=False)

			if translate:
				to_r.index = row2id_batch( to_r.index.tolist(), host, port, db, return_field = "egrin2_row_name" )

		if y_type == "columns":
			cols = pd.Series( list( itertools.chain( *query["columns"].tolist() ) ) ).value_counts().to_frame( "counts" )
			
			# filter out rows that aren't in the database - i.e. not annotated in MicrobesOnline
			in_db = pd.DataFrame( list( client[db].col_info.find( { }, { "_id" : 0, "col_id": 1 } ) ) ).col_id.tolist()
			common_cols = list( set( cols.index ).intersection( set( in_db )  ))
			cols = cols.loc[common_cols]
			
			# find all bicluster counts
			all_counts = pd.DataFrame( list( client[db].bicluster_info.find( { }, { "columns": 1 } ) ) )
			all_counts = pd.Series( list( itertools.chain( *all_counts["columns"].tolist() ) ) ).value_counts().to_frame( "all_counts" )
			# combine two data frames
			to_r = cols.join(all_counts).sort("counts",ascending=False)

			if translate:
				to_r.index = col2id_batch( to_r.index.tolist(), host, port, db, return_field = "egrin2_col_name" )


		if y_type == "motif.gre_id":
			gres = pd.Series( list( itertools.chain( *[ i.values() for i in list( itertools.chain( *query.motif.tolist() ) ) ] ) ) )
			gres = gres.value_counts().to_frame( "counts" )

			# find all bicluster counts
			all_counts = pd.DataFrame( list( client[db].bicluster_info.find( { }, { "motif.gre_id": 1 } ) ) )
			all_counts = pd.Series( list( itertools.chain( *[ i.values() for i in list( itertools.chain( *all_counts.motif.tolist() ) ) ] ) ) ).value_counts().to_frame( "all_counts" )
			# combine two data frames
			to_r = gres.join(all_counts).sort("counts",ascending=False)

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