#!/usr/bin/env python

"""Resample ratios. Add them to col_resample collection in an existing egrin2 MongoDB"""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import random
from multiprocessing import Pool

# $ hg clone ssh://hg@bitbucket.org/djcbeach/monary ./monary
# $ cd ./monary && python setup.py install
from monary import Monary

from pymongo import MongoClient
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

def multiSample(number):
	"""Returns how many documents in the collection"""
	colResampleInd( rows = rows, cols = [server_number])
 
def rsd( vals ):
	return abs( np.std( vals ) / np.mean( vals ) )

def resample( row_vals, n_rows ):
		return rsd( random.sample( row_vals, n_rows ) )

def choose_n( col, vals, n, add, client, db, n_rows, n_resamples ):
	normalized = vals.loc[:,"normalized_expression"].copy()
	normalized.sort()
	standardized = vals.loc[:,"standardized_expression"].copy()
	standardized.sort()

	if add:
		d = {
			"n_rows": n_rows,
			"col_id": col,
			"resamples": n_resamples,
			"lowest_normalized": normalized.iloc[0:n].tolist(),
	 		"lowest_standardized": standardized.iloc[0:n].tolist()
		}
		client[ db ][ "col_resample" ].insert( d )
	else:
		#update
		resamples = n_resamples + old_records[ col ][ "resamples" ]
		n2keep = round( resamples*keepP )
		ran = normalized + old_records[ col ][ "lowest_normalized" ]
		ran.sort()
		ran = ran.iloc[ 0: int( n2keep ) ].tolist()
		ras = standardized + old_records[ col ][ "lowest_standardized" ]
		ras.sort()
		ras = ras.iloc[ 0: int( n2keep ) ].tolist()
		client[ db ][ "col_resample" ].update( { "n_rows": n_rows, "col_id": col }, { "$set": { "resamples": resamples, "lowest_normalized": ran, "lowest_standardized": ras } } )

def colResampleInd( host, n_rows, cols, n_resamples = 20000, keepP = 0.1, port = 27017, db = "egrin2_db" ):
	"""Resample gene expression for a given number of genes in a particular condition using RSD"""

	print n_rows

	# make connection
	client = MongoClient( host = host, port=port )

	# see if a record for this gene set size exists 
	old_records = { i[ "col_id" ] : i for i in client[ db ][ "col_resample" ].find( { "n_rows": n_rows, "col_id": { "$in": cols } } ) }

	if  old_records is not None:
		toUpdate = [ int( i[ "col_id" ] ) for i in old_records.values() if int( i[ "resamples" ] ) < n_resamples ]
		
	toAdd = [ i for i in cols if i not in old_records.keys() ]

	if len( toAdd ) == 0 and len( toUpdate ) == 0:
		print "Nothing to add"
		client.close()
		return None

	n2keep = int( round( n_resamples*keepP ) )

	# toAdd
	if len( toAdd ) > 0:
		print "Adding new document to MongoDB"
		with Monary( 'mongodb://'+host+':'+str(port)+'/' ) as mon:
			mon.connect()
			columns = ["col_id", "normalized_expression", "standardized_expression"]
			arrays = mon.query(
				"egrin2_db",                         # database name
				"gene_expression",                   # collection name
				{ "col_id": { "$in": toAdd } },                             # query spec
				columns, # field names (in Mongo record)
				["int64", "float64", "float64"] 
			)
			df = np.matrix(arrays).transpose()
			df = pd.DataFrame(df, columns=columns)

		df = df.groupby("col_id")
		df_rsd = pd.concat( [ df.aggregate( resample, n_rows ) for i in range( 0, n_resamples ) ] )
		df_rsd = df_rsd.groupby( df_rsd.index )
		tmp = [ choose_n( int( i ), df_rsd.get_group( i ), n2keep, True, client, db, n_rows, n_resamples ) for i in df_rsd.groups.keys() ]

	# toUpdate
	if len( toUpdate ) > 0:
		print "Updating existing MongoDB documents"
		with Monary( 'mongodb://'+host+':'+str(port)+'/' ) as mon:
			mon.connect()
			columns = ["col_id", "normalized_expression", "standardized_expression"]
			arrays = mon.query(
				"egrin2_db",                         # database name
				"gene_expression",                   # collection name
				{ "col_id": { "$in": toUpdate } },                             # query spec
				columns, # field names (in Mongo record)
				["int64", "float64", "float64"] 
			)
			df = np.matrix(arrays).transpose()
			df = pd.DataFrame(df, columns=columns)

		df = df.groupby("col_id")
		resamples = n_resamples - np.min( [ i[ "resamples" ] for i in old_records.values( ) ] )
		if resamples > 0:
			df_rsd = pd.concat( [ df.aggregate( resample, n_rows ) for i in range( 0, resamples ) ] )
			df_rsd = df_rsd.groupby( df_rsd.index )
			tmp = [ choose_n( int( i ), df_rsd.get_group( i ), n2keep, False, client, db, n_rows, resamples ) for i in df_rsd.groups.keys() ]
	client.close()

	return None

if __name__ == '__main__':

	host = "localhost"
	port = 27017
	db = "egrin2_db"

	client = MongoClient( host = host, port=port )
	
	cols = range( 0,client[ db ][ "col_info" ].count( ) )

	tmp = Parallel(n_jobs=4)( delayed( colResampleInd )( "localhost", i, cols, n_resamples = 10) for i in range( 3,10 ) )

	print "Done"	