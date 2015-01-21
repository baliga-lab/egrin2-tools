#!/usr/bin/env python

"""Tools for querying EGRIN2.0 MongoDB."""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks"]
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

def rsd( vals ):
	return abs( np.std( vals ) / np.mean( vals ) )

def check_colResamples( col, n_rows, n_resamples, host = "localhost", port = 27017, db = "egrin2_db",  ):
	# connect to db
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' ) 
	if client[db].col_resample.find_one( { "n_rows": n_rows, "col_id": col, "resamples": { "$gte": n_resamples } } ) is None:
		client.close()
		return col
	client.close()

def row2id( row, host = "localhost", port = 27017, db = "egrin2_db",  verbose = False ):
	"""Check name format of rows. If necessary, translate."""
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = list( client[db].row_info.find( { "$or": [ { "row_id": row }, { "egrin2_row_name": row }, { "GI": row }, { "accession": row }, { "name": row }, { "sysName": row } ] } ) )
	client.close()
	if len( query ) == 1: 
		row = query[ 0 ][ "row_id" ]
		return row
	elif len( query ) > 0:
		print "ERROR: Multiple genes match the row name: %s" % row
		if verbose:
			print query
		return None
	else:
		print "ERROR: Cannot identify row name: %s" % row
		return None
	

def col2id( col, host = "localhost", port = 27017, db = "egrin2_db",  verbose = False ):
	"""Check name format of rows. If necessary, translate."""
	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	query = list( client[db].col_info.find( { "$or": [ { "col_id": col }, { "egrin2_col_name": col } ] } ) )
	client.close()
	if len( query ) == 1: 
		col = query[ 0 ][ "col_id" ]
		return col
	elif len( query ) > 0:
		print "ERROR: Multiple genes match the row name: %s" % col
		if verbose:
			print query
		return None
	else:
		print "ERROR: Cannot identify row name: %s" % col
		return None

def col2name( col, host = "localhost", port = 27017, db = "egrin2_db",  verbose = False ):
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

def rowsColPval( rows = None, cols = None, host = "localhost", port = 27017, db = "egrin2_db", standardized = None, sig_cutoff = None, sort = True ):

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

	rows = [ row2id( i, host, port, db ) for i in rows ]
	rows = [ i for i in rows if i is not None]
	if rows is None:
		print "Please provide an appropriately named array of rows"
		return None 
	
	cols = [ col2id( i, host, port, db ) for i in cols ]
	cols = [ i for i in cols if i is not None]
	if cols is None:
		print "Please provide an appropriately named array of cols"
		return None

	if standardized is None:
		# compute RSD on standardized gene expression by default
		# other option 'False' for normalized (not standardized expression)
		standardized = True

	if sig_cutoff is None:
		# return only cols equal below sig_cutoff
		sig_cutoff = 0.05

	client = MongoClient( 'mongodb://'+host+':'+str(port)+'/' )
	exp_df = pd.DataFrame( list( client[ db ].gene_expression.find( { "col_id": { "$in" : cols }, "row_id": { "$in" : rows } }, { "_id":0, "col_id":1, "normalized_expression":1, "standardized_expression":1 } ) ) )
	random_rsd = pd.DataFrame( list( client[ db ].col_resample.find( { "n_rows": len( rows ), "col_id": { "$in" : cols } }, { "_id":0 } ) ) )
	random_rsd.index = random_rsd["col_id"]

	if random_rsd is None:
		print "Could not find resample DB entry for %i rows" % ( len(rows) )
		return None

	exp_df_rsd = exp_df.groupby("col_id").aggregate(rsd)

	if standardized:
		exp_df_rsd = exp_df_rsd.loc[ :, "standardized_expression"]
		resamples = random_rsd.loc[ :, "resamples"].to_dict()
		random_rsd = random_rsd.loc[ :,"lowest_standardized" ].to_dict()
	else:
		exp_df_rsd = exp_df_rsd.loc[ :, "normalized_expression"]
		resamples = random_rsd.loc[ :, "resamples"].to_dict()
		random_rsd = random_rsd.loc[ :,"lowest_normalized" ].to_dict()
		

	pvals = exp_df_rsd.groupby(level=0).aggregate(empirical_pval, random_rsd, resamples )
	pvals.columns = ["pval"]
	pvals.index = [ col2name( i ) for i in pvals.index.values]

	if sig_cutoff is not None:
		pvals = pvals[ pvals <= sig_cutoff ]

	if sort:
		pvals.sort()

	pvals = pvals.to_frame()
	pvals.columns = [ "pval" ] 

	return pvals
	
def colResampleGroup( self, rows = None, cols = None, n_resamples = None, host = "localhost", port = 27017, db = "egrin2_db", sig_cutoff = None, standardized = None ):
	"""Resample gene expression for a given set of genes in any number of conditions. Should be used instead of colReampleInd (it calls that function)"""

	if n_resamples is None:
		n_resamples = 20000

	if standardized is None:
		# compute RSD on standardized gene expression by default
		# other option 'False' for normalized (not standardized expression)
		standardized = True

	if sig_cutoff is None:
		# return only cols equal below sig_cutoff
		sig_cutoff = 0.05

	rows = [ row2id( i, host, port, db ) for i in rows ]
	rows = [ i for i in rows if i is not None]
	if rows is None:
		print "Please provide an appropriately named array of rows"
		return None 
	
	cols = [ col2id( i, host, port, db ) for i in cols ]
	cols = [ i for i in cols if i is not None]
	if cols is None:
		print "Please provide an appropriately named array of cols"
		return None
	
	# Determine what/how many resamples need to be added to db
	toAdd = [ check_colResamples( i, len( rows ), n_resamples, host, port , db ) for i in cols]
	toAdd = [ i for i in toAdd if i is not None]

	count = 1
	if len( toAdd) > 0:
		print "I need to perform %i random resample(s) of size %i to compute pvals. Please be patient. This may take a while..." % ( len(toAdd), n_resamples )
		for i in toAdd:		
			currentEntry = self.db.col_resample.find_one( { "n_rows": len( rows ), "col_id": i } )
			if currentEntry is not None:
				# only add enough resamples to reach requested value
				self.colResampleInd( len( rows ), i, n_resamples - currentEntry[ "resamples" ], .1 )
			else:
				self.colResampleInd( len( rows ), i, n_resamples, .1 )
			if round( ( float( count) / len(toAdd) ) * 100, 2 ) % 10  == 0:
				print "%d percent" % ( round( ( float( count ) / len( toAdd ) ) * 100, 2 ) )
			count = count + 1
		print "Done adding random resamples."

	print "Calculating pvals"

	pvals = self.rowsColPval( rows, cols, standardized, sig_cutoff )

	return pvals

if __name__ == '__main__':

	host = "primordial"
	port = 27017
	db = "egrin2_db"