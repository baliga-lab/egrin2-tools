def checkDBresamples( self, col, n_rows, n_resamples ):
		if self.db.col_resample.find_one( { "n_rows": n_rows, "col_id": col, "resamples": { "$gte": n_resamples } } ) is None:
			return col

def checkRow( self, row ):
	"""Check name format of rows. If necessary, translate."""
	if row in self.row2id.keys():
		return self.row2id[ row ] 
	elif int( row ) in self.id2row.keys():
		return int( row )
	else:
		print "ERROR: Cannot identify row name: %s" % row
		return None

def checkCol( self, col ):
	"""Check name format of rows. If necessary, translate."""
	if col in self.col2id.keys():
		return self.col2id[ col ] 
	elif int( col ) in self.id2col.keys():
		return int( col )
	else:
		print "ERROR: Cannot identify col name: %s" % col
		return None

def rowsColPval( self, rows = None, cols = None, standardized = None, sig_cutoff = None ):

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

	rows = [ self.checkRow( i ) for i in rows ]
	rows = [ i for i in rows if i is not None]
	if rows is None:
		print "Please provide an appropriately named array of rows"
		return None 
	
	cols = [ self.checkCol( i ) for i in cols ]
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

	exp_df = pd.DataFrame( )

	# compute actual RSD
	if standardized:
		 exp_df =  exp_df.append( [ pd.Series( [ i[ "col_id" ], i[ "standardized_expression" ] ], index = [ "col_id", "value" ] ) for i in self.db.gene_expression.find( { "col_id": { "$in" : cols }, "row_id": { "$in" : rows } } ) if isinstance( i[ "standardized_expression" ], float ) ], ignore_index = True )
		 # get random resamples for this col
		 random_rsd = { i[ "col_id" ]: i[ "lowest_standardized" ] for i in self.db.col_resample.find( { "n_rows": len( rows ), "col_id": { "$in" : cols } } ) }
		 resamples = { i[ "col_id" ]: i[ "resamples" ] for i in self.db.col_resample.find( { "n_rows": len( rows ), "col_id": { "$in" : cols } } ) }
	else:
		exp_df =  exp_df.append( [ pd.Series( [ i[ "col_id" ], i[ "normalized_expression" ] ], index = [ "col_id", "value" ] ) for i in self.db.gene_expression.find( { "col_id": { "$in" : cols }, "row_id": { "$in" : rows } } ) if isinstance( i[ "normalized_expression" ], float ) ], ignore_index = True )
		# get random resamples for this col
		random_rsd = { i[ "col_id" ]: i[ "lowest_normalized" ] for i in self.db.col_resample.find( { "n_rows": len( rows ), "col_id": { "$in" : cols } } ) }
		resamples = { i[ "col_id" ]: i[ "resamples" ] for i in self.db.col_resample.find( { "n_rows": len( rows ), "col_id": { "$in" : cols } } ) }

	if random_rsd is None:
		print "Could not find resample DB entry for %i rows" % ( len(rows) )
		return None

	exp_df_rsd = exp_df.groupby("col_id").aggregate(self.rsd)

	pvals = exp_df_rsd.groupby(level=0).aggregate(empirical_pval, random_rsd, resamples )
	pvals.columns = ["pval"]
	pvals.index = [ self.id2col[ i ] for i in pvals.index.values]

	if sig_cutoff is not None:
		pvals = pvals[ pvals.loc[ :,"pval" ] <= sig_cutoff ]

	return pvals
	
def colResampleGroup( self, rows = None, cols = None, n_resamples = None, sig_cutoff = None, standardized = None ):
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

	rows = [ self.checkRow( i ) for i in rows ]
	rows = [ i for i in rows if i is not None]
	if rows is None:
		print "Please provide an appropriately named array of rows"
		return None 
	
	cols = [ self.checkCol( i ) for i in cols ]
	cols = [ i for i in cols if i is not None]
	if cols is None:
		print "Please provide an appropriately named array of cols"
		return None
	
	# Determine what/how many resamples need to be added to db
	toAdd = [ self.checkDBresamples( i, len( rows ), n_resamples ) for i in cols]
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