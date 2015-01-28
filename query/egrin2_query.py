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

def col2name( col, host, port, db,  verbose = False ):
	"""Check name format of rows. If necessary, translate."""
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = list( client[db].col_info.find( { "$or": [ { "col_id": col }, { "egrin2_col_name": col } ] } ) )
	client.close()
	if len( query ) == 1: 
		col = query[ 0 ][ "egrin2_col_name" ]
		return col
	elif len( query ) > 0:
		print "ERROR: Multiple genes match the row name: %s" % col
		if verbose:
			print query
		return None
	else:
		print "ERROR: Cannot identify row name: %s" % col
		return None

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

def x2bicluster( x = [ 0,1 ], x_type = None, logic = "and", count = True, host = "localhost", port = 27017, db = "",  verbose = False, return_field = [ "cluster" ] ):
	"""
	Determine how often 'x' occurs in biclusters. Usually just retrieve the counts. Retrieve additional bicluster info by setting count to False

	Available x_type(s):

	'rows' or 'genes': search for row (genes) in biclusters. x should be a list of rows, eg ["carA","carB"] or [275, 276]
	'columns' or 'conditions': search for columns (conditions) in biclusters. x should be a list of columns, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
	'gre': search for GREs in biclusters. x should be a list of GRE IDs, eg [4, 19]
	''

	"""

	if x_type is None:
		print "Please supply a type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres' "

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	if x_type == "rows" or "row" or "gene" or "genes":
		x_type = "rows"
		x_o = x
		x = [ row2id( i, host, port, db ) for i in x ]
		x = [ i for i in x if i is not None]
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate row names: %s" % x_o
			return None

	if x_type == "columns" or "column" or "condition" or "conditions":
		x_type = "columns"
		x_o = x
		x = [ col2id( i, host, port, db ) for i in x ]
		x = [ i for i in x if i is not None]
		x = list( set( x ) )
		if len( x ) == 0:
			print "Cannot translate col names: %s" % x_o
			return None  

	if x_type == "motif" or "gre" or "motc" or "motif.gre" or "motfs" or "gres" or "motcs":
		x_type = "motif.gre_id"

	if x_type in [ "rows", "columns", "motif.gre_id" ] and logic in [ "and","or","nor" ]:
		q = { "$"+logic: [ { x_type : i } for i in x ] }
		if count:
			query =client[db].bicluster_info.find( q ).count()
			return query
		else:
			query = pd.DataFrame( list( client[db].bicluster_info.find( q ) ) )
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
		print "Could not find any biclusters matching your criteria"
		return None

def bicluster2x( bicluster_ids = [ ], x_type = None, host = "localhost", port = 27017, db = "",  verbose = False ):
	"""
	Determine how often 'x' occurs in biclusters. Usually just retrieve the counts. 

	Available x_type(s):

	'rows' or 'genes': search for row (genes) in biclusters. x should be a list of rows, eg ["carA","carB"] or [275, 276]
	'columns' or 'conditions': search for columns (conditions) in biclusters. x should be a list of columns, eg ["dinI_U_N0025", "dinP_U_N0025"] or [0,1]
	'gre': search for GREs in biclusters. x should be a list of GRE IDs, eg [4, 19]
	''

	"""
	if x_type is None:
		print "Please supply a type for your query. Types include: 'rows' (genes), 'columns' (conditions), 'gres' "

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )

	if x_type == "rows" or "row" or "gene" or "genes":
		x_type = "rows"

	if x_type == "columns" or "column" or "condition" or "conditions":
		x_type = "columns"

	if x_type == "motif" or "gre" or "motc" or "motif.gre" or "motfs" or "gres" or "motcs":
		x_type = "motif.gre_id"

	if x_type in [ "rows", "columns", "motif.gre_id" ] and logic in [ "and","or","nor" ]:
		q = { "_id": { "$in": bicluster_ids } }
		o = { x_type: 1 }
		query = pd.DataFrame( list( client[db].bicluster_info.find( q, o ) ) )
	else:
		print "I don't recognize the logic you are trying to use. 'logic' must be 'and', 'or', or 'nor'."
	
	client.close()

	if query.shape[0] > 0: 
		if x_type == "rows":
			in_db = pd.DataFrame( list( client[db].row_info.find( { }, { "_id" : 0, "row_id": 1 } ) ) ).row_id.tolist()
			rows = pd.Series( list( itertools.chain( *query.rows.tolist() ) ) ).value_counts().to_frame( "counts" )
			# filter out rows that aren't in the database - i.e. not annotated in MicrobesOnline
			common_rows = list(set(rows.index).intersection(set(in_db)))
			rows = rows.loc[common_rows]

			rows[ "pval" ] = 
			t_rows = [ row2id( i, host, port, db,  verbose = False, return_field = "egrin2_row_name" ) for i in rows.index ]
			rows.index = 
		
		if return_field == "all":
			return query
		else:
			try:
				return query.loc[ :, return_field ]
			except Exception:
				return query
	else:
		print "Could not find any biclusters matching your criteria"
		return None

	return None

def agglom( x, x_type, out_type, path ):
	return None



	
	
	
	

if __name__ == '__main__':

	host = "primordial"
	port = 27017
	db = "egrin2_db"