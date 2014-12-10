#!/usr/bin/env python

"""Choose ensemble conditions"""

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
import random

import pdb

import numpy as np
import pandas as pd

import pymongo
from pymongo import MongoClient
import gridfs
from Bio import SeqIO

# ratios = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141130.MTB.all.ratios.csv"
# blocks = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141202.MTB.EGRIN2.blocks.csv"
# exclusion = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141202.MTB.EGRIN2.exclusion.blocks.csv"
# inclusion = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141202.MTB.EGRIN2.inclusion.blocks.csv"
# tfblocks = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/TF.blocks.csv"

# ensemblePicker( ratios, blocks, exclusion, inclusion, tfblocks = None, nruns=500, exclusion_percentage=0.25 )
random.normalvariate


# tmp.blocks.block.value_counts()

class ensemblePicker:
    	"""Pick conditions for an ensemble run, biasing towards inclusion of blocks of conditions in the inclusion blocks while making sure that conditions in exclusion blocks are excluded together in at least some percentage of runs"""
	def __init__( self, ratios, blocks, exclusion, inclusion, avg_col_size = None, sd_col_size = None, inclusion_weight = None, tfblocks = None, nruns = None, exclusion_percentage = None, n_rand_exclusion = None ):
		
		if nruns == None:
			self.nruns = 100
		
		if exclusion_percentage == None:
			self.exclusion_percentage = 0.25

		if inclusion_weight == None:
			self.inclusion_weight = 2

		if avg_col_size == None:
			self.avg_col_size = 250

		if sd_col_size == None:
			self.sd_col_size = 100

		if n_rand_exclusion == None:
			self.n_rand_exclusion = 3

		self.ratios = pd.read_csv( ratios, index_col=0, sep="," )
		self.blocks2col = pd.read_csv( blocks, sep="," ).icol([0,1])
		self.blocks =  pd.DataFrame( zip( self.blocks2col.block.value_counts().keys(),self.blocks2col.block.value_counts() ), columns = ["blocks", "block.sample.num"]) 
		self.blocks["p"] = 1
		self.exclusion = pd.read_csv( exclusion, sep="," )
		self.exclusion["p"] = 1
		self.inclusion = pd.read_csv( inclusion, sep="," )
		self.inclusion["p"] = 1
		self.inclusion_matrix = pd.DataFrame( 1, index = np.unique(self.blocks2col.block) , columns = np.unique(self.blocks2col.block) )

	def chooseRandomBlocks( self, blocks, exclusion, n_rand_exclusion ):
		# maximum block size
		max_b_size = round( .10  * self.ratios.shape[1] )
		names = []
		n_cols = []
		for i in range( 0, n_rand_exclusion ):
			n_c = 0
			c = []
			while n_c < max_b_size:
				# choose block
				j = random.choice( range( 0, self.blocks.shape[0] ) )
				c.append( self.blocks.iloc[ j ][ "blocks" ] ) 
				n_c = n_c + self.blocks.iloc[ j ][ "block.sample.num" ]
			names.append( ":::".join(c) )
			n_cols.append( n_c )
		tmp_df = pd.DataFrame( zip( names, n_cols ), columns = [ "exclusion.blocks", "block.sample.num" ] )
		tmp_df[ "p" ] = 1
		df = pd.concat( [ exclusion,tmp_df ], ignore_index=True )
		return df


		
	def weightedReservoirSample( self, n_samples, sample_names, weights ):
		"""Choose a block loosely based on weighted reservoir sampling - the easy way"""
		k = [np.power( random.random(), 1.0/weights[i] ) for i in sample_names ]
		df = pd.DataFrame( zip( sample_names, k ), index = sample_names, columns = [ "sample_name", "key"] ).sort("key",ascending = False )
		k.sort( reverse = True )
		return( df.iloc[ 0 : n_samples ].index.values )

	def updateWeights(self, ):
		"""Update probability of picking conditions based on its current representation in the ensemble"""

	def pickCols_single( self ):
		"""Pick columns for a cMonkey run using predefined blocks and based on their current representation in the ensemble"""
		#
		tmp = []

	def pickCols_all( self ):
		tmp = []

