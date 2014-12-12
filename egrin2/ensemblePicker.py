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
import itertools

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

# ensemblePicker( ratios, blocks, exclusion, inclusion, tfblocks = None, nruns=100, exclusion_percentage=0.25 )

# tmp.blocks.block.value_counts()

class ensemblePicker:
    	"""Pick conditions for an ensemble run, biasing towards inclusion of blocks of conditions in the inclusion blocks while making sure that conditions in exclusion blocks are excluded together in at least some percentage of runs"""
	def __init__( self, ratios, blocks, exclusion, inclusion, avg_col_size = None, sd_col_size = None, inclusion_weight = None, tfblocks = None, nruns = None, exclusion_percentage = None, n_rand_exclusion = None ):
		
		if nruns == None:
			self.nruns = 100
		else:
			self.nruns = nruns
		
		if exclusion_percentage == None:
			self.exclusion_percentage = 25
		else:
			self.exclusion_percentage = exclusion_percentage

		if inclusion_weight == None:
			self.inclusion_weight = 2
		else:
			self.inclusion_weight = inclusion_weight

		if avg_col_size == None:
			self.avg_col_size = 250
		else:
			self.avg_col_size = avg_col_size

		if sd_col_size == None:
			self.sd_col_size = 100
		else:
			self.sd_col_size = sd_col_size

		def strip(text):
		    try:
		        return text.strip()
		    except AttributeError:
		        return text

		self.ratios = pd.read_csv( ratios, index_col=0, sep="," )
		self.blocks2col = pd.read_csv( blocks, sep=",", names=[ "sample", "block" ], converters = {'sample' : strip,
                                    'block' : strip,
                                    } ).icol( [0,1] )
		self.blocks =  pd.DataFrame( zip( self.blocks2col.block.value_counts().keys(), self.blocks2col.block.value_counts() ), columns = ["block","block.sample.num"] ) 
		self.blocks["p"] = 1
		self.blocks["r_in"] = 0
		self.blocks = self.blocks.set_index(keys="block")
		self.blocks = self.blocks[self.blocks.index != "block"]
		self.blocks2col = self.blocks2col.set_index(keys="block")
		self.exclusion = pd.read_csv( exclusion, sep=",", converters = { "exclusion.blocks" : strip }, index_col = 0 )
		self.exclusion["p"] = 1
		self.exclusion["r_out"] = 0
		self.inclusion = pd.read_csv( inclusion, sep=",", converters = { "inclusion.blocks" : strip }, index_col = 0 )
		self.inclusion["p"] = 1
		self.inclusion["r_in"] = 0
		self.inclusion_matrix = pd.DataFrame( 0, index = self.blocks.index , columns = self.blocks.index )

		# fill out inclusion matrix
		for i in self.inclusion.index:
			blocks = i.split( ":::" )
			for i, j in itertools.combinations( blocks, 2 ):
				i = i.strip()
				j = j.strip()
				self.inclusion_matrix.loc[ i, j ] = self.inclusion_weight
				self.inclusion_matrix.loc[ j, i ] = self.inclusion_weight

		if n_rand_exclusion == None:
			self.n_rand_exclusion = round( self.nruns / self.exclusion_percentage ) - self.exclusion.shape[0]
		else:
			self.n_rand_exclusion = n_rand_exclusion

		if self.n_rand_exclusion < 0:
			print "You have defined more exclusion blocks than can be excluded given the number of runs at the requested exclusion rate. Maximum exclusion rate for %s runs is %s percent!" % (self.nruns, self.nruns / self.exclusion.shape[0])

		self.run_composition = {}	

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
		tmp_df[ "r_out" ] = 0
		df = pd.concat( [ exclusion,tmp_df ], ignore_index=True )
		return df
		
	def weightedReservoirSample( self, n_samples, sample_names, weights ):
		"""Choose a block loosely based on weighted reservoir sampling - the easy way"""
		k = [np.power( random.random(), 1.0/weights[i] ) for i in sample_names ]
		df = pd.DataFrame( zip( sample_names, k ), index = sample_names, columns = [ "sample_name", "key"] ).sort("key",ascending = False )
		k.sort( reverse = True )
		return( df.iloc[ 0 : n_samples ].index.values )

	def updateWeights( self ):
		"""Update probability of picking blocks based on their current representation in the ensemble"""
		self.blocks[ "p" ] = ( 1+10e-10 ) - self.blocks[ "r_in" ]
		self.inclusion[ "p" ] = ( 1+10e-10 ) - self.inclusion[ "r_in" ]
		self.exclusion[ "p" ] = np.power( ( 1+10e-10 ) - self.exclusion[ "r_out" ], abs( self.exclusion_percentage/100.0 - self.exclusion[ "r_out" ] ) / 1 )

	def combinedWeights( self, excluded_diff, blocks ):
		tmp_m = self.inclusion_matrix.loc[ excluded_diff, excluded_diff ]
		tmp_m = tmp_m.loc[ set( blocks ), ].sum(0)
		tmp_m[ tmp_m==0] = 1
		return tmp_m

	def blocks2cols(self, blocks, excluded ):
		"""Get col ids from blocks. Since a single condition can be assigned to multiple blocks, we need to make sure that conditions in a current exlcusion block did not 'sneak in' to the run"""
		block_cols = self.blocks2col.loc[ blocks.split(":::") ] 
		excluded_cols = self.blocks2col.loc[ excluded.split(":::") ] 
		return set( block_cols.sample ).difference( set( excluded_cols.sample ) )

	def pickCols_single( self, n ):
		"""Pick columns for a cMonkey run using predefined blocks and based on their current representation in the ensemble"""
		# what is max conditions for this run?
		n_cols = self.nruns+1
		while n_cols > self.nruns:
			n_cols = round( self.avg_col_size + ( self.avg_col_size/4 )*random.gammavariate(1,2) )
		
		# first choose excluded blocks
		excluded = self.weightedReservoirSample( 1, self.exclusion.index, self.exclusion["p"] )[0]
		# update r_out
		self.exclusion.loc[ excluded, "r_out" ] = 1.0/self.nruns + self.exclusion.loc[ excluded, "r_out" ] 

		excluded_diff = set( self.blocks.index ).difference( set( excluded.split( ":::" ) ) )

		# pick a single from blocks remaining
		blocks = []
		cols = []
		while len(cols) < n_cols:
			print blocks
			if len(cols) == 0:
				# this is the first condition
				weights = self.blocks.loc[ excluded_diff, "p" ]
				block_one = self.weightedReservoirSample( 1, excluded_diff, weights )[0] 
				col_one = self.blocks2cols( block_one, excluded )
				blocks.append( block_one )
				[cols.append( i ) for i in col_one]
				excluded_diff = excluded_diff.difference( set( [ block_one ] ) )
			else:
				weights = self.blocks.loc[ excluded_diff, "p" ] * self.combinedWeights( excluded_diff, blocks )
				block_one = self.weightedReservoirSample( 1, excluded_diff, weights )[0] 
				col_one = self.blocks2cols( block_one, excluded )
				blocks.append( block_one )
				[cols.append( i ) for i in col_one]
				excluded_diff = excluded_diff.difference( set( [ block_one ] ) )
				print len( excluded_diff )

		self.run_composition[n] = {
		"blocks" : blocks,
		"cols": cols,
		"excluded": excluded
		}
		return None


		

	def pickCols_all( self ):
		# get random blocks
		self.exclusion = chooseRandomBlocks( self.blocks, self.exclusion, self.n_rand_exclusion )
		tmp = []

