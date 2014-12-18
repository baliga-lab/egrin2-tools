#!/usr/bin/env python

"""Select ensemble conditions given defined blocks of conditions to include together (inclusion blocks) or exclude together (exclusion blocks)"""

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

# eg. how to run

# from egrin2.ensemblePicker import *

# ratios = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141130.MTB.all.ratios.csv"
# blocks = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141202.MTB.EGRIN2.blocks.csv"
# exclusion = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141202.MTB.EGRIN2.exclusion.blocks.csv"
# inclusion = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/20141202.MTB.EGRIN2.inclusion.blocks.csv"
# tfblocks = "/Users/abrooks/Dropbox/MTB.EGRIN2.0.files/TF.blocks.csv"

# tmp = ensemblePicker( ratios, blocks, exclusion, inclusion, tfblocks = None, nruns=500, exclusion_percentage=25 )
# tmp.pickCols_all()

class ensemblePicker:
    	"""Pick conditions for an ensemble run, biasing towards inclusion of blocks of conditions in the inclusion blocks while making sure that conditions in exclusion blocks are excluded together in at least some percentage of runs"""
	def __init__( self, ratios, blocks, exclusion, inclusion, avg_col_size = None, sd_col_size = None, inclusion_weight = None, tfblocks = None, nruns = None, exclusion_percentage = None, n_rand_exclusion = None, report_file = None, random_blocks = None ):
		
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

		if report_file == None:
			self.report_file = None

		if random_blocks == None:
			self.random_blocks = True

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
		self.blocks2col["r_in"] = 0

		self.exclusion = pd.read_csv( exclusion, sep=",", converters = { "exclusion.blocks" : strip }, index_col = 0 )
		self.exclusion["p"] = 1
		self.exclusion["r_out"] = 0
		self.inclusion = pd.read_csv( inclusion, sep=",", converters = { "inclusion.blocks" : strip }, index_col = 0 )
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
			self.n_rand_exclusion = int( 1/ ( float( self.exclusion_percentage )/100.0 ) ) - self.exclusion.shape[0]
			if self.n_rand_exclusion < 0:
				self.n_rand_exclusion = 0
		else:
			self.n_rand_exclusion = n_rand_exclusion

		if self.n_rand_exclusion < 0:
			print "You have defined more exclusion blocks than can be excluded given the number of runs at the requested exclusion rate. Maximum exclusion rate for %s runs is %s percent!" % (self.nruns, self.nruns / self.exclusion.shape[0])

		self.run_composition = {}	

	def chooseRandomBlocks( self ):
		# maximum block size
		max_b_size = round( .10  * self.ratios.shape[1] )
		names = []
		n_cols = []
		for i in range( 0, int( self.n_rand_exclusion ) ):
			n_c = 0
			c = []
			while n_c < max_b_size:
				# choose block
				j = random.choice( range( 0, self.blocks.shape[0] ) )
				c.append( self.blocks.iloc[ j ].name ) 
				n_c = n_c + self.blocks.iloc[ j ][ "block.sample.num" ]
			names.append( ":::".join(c) )
			n_cols.append( n_c )
		tmp_df = pd.DataFrame( n_cols, index = names, columns = [ "block.sample.num" ] )
		tmp_df[ "p" ] = 1
		tmp_df[ "r_out" ] = 0
		self.exclusion = pd.concat( [ self.exclusion, tmp_df ] )
		return None
		
	def weightedReservoirSample( self, n_samples, sample_df ):
		"""Choose a block loosely based on weighted reservoir sampling - the easy way"""
		k = pd.DataFrame([np.power( random.random(), 1.0/sample_df.loc[i] ) for i in sample_df.index ], index = sample_df.index, columns = ["p"]).sort("p",ascending = False)
		return( k.iloc[ 0 : n_samples ].index.values )

	def updateWeights( self, blocks ):
		"""Update probability of picking blocks based on their current representation in the ensemble"""
		#print "updating weights"
		self.blocks.loc[ blocks, "r_in" ] = self.blocks.loc[ blocks, "r_in" ] + 1.0/self.nruns
		self.blocks[ "p" ] = ( 1+10e-10 ) - self.blocks[ "r_in" ]

	def combinedWeights( self, excluded_diff, blocks ):
		tmp_m = self.inclusion_matrix.loc[ blocks, excluded_diff ].sum(0)
		tmp_m[ tmp_m==0] = 1
		return tmp_m

	def translate_blocks(self, blocks, excluded ):
		"""Get col ids from blocks. Since a single condition can be assigned to multiple blocks, we need to make sure that conditions in a current exclusion block did not 'sneak in' to the run"""
		block_cols = self.blocks2col.loc[ blocks.split(":::") ] 
		excluded_cols = self.blocks2col.loc[ excluded.split(":::") ] 
		return set( block_cols.sample ).difference( set( excluded_cols.sample ) )

	def count_cols(self, col):
		self.blocks2col.loc[ col, "r_in"] = self.blocks2col.loc[ col, "r_in"] + 1.0/self.nruns

	def findBlockInclusionFreq( self ):
		blocks = self.inclusion.index
		block_stats = pd.DataFrame(0, index = blocks, columns = [ "freq_single", "freq_coinclusion", "max_coinclusion_all", "max_coinclusion_sub" ] )
		for block in blocks:
			block_ind = block.split(":::") 
			co_freq = pd.DataFrame(0, index = block_ind, columns = self.run_composition.keys() )
			for i in self.run_composition.iterkeys():
				for j in block_ind:
					if j in self.run_composition[i]["blocks"]:
						co_freq.loc[ j, i ] = co_freq.loc[ j, i ] + 1
			block_stats.loc[ block, "freq_single" ] = round( sum( co_freq.sum( 0 ) > 0 ) / float( co_freq.shape[1] ), 2 )
			block_stats.loc[ block, "freq_coinclusion" ] = round( sum( co_freq.sum(1) )/float( ( co_freq.shape[0] * sum(co_freq.sum( 0 ) > 0) ) ), 2)
			block_stats.loc[ block, "max_coinclusion_all" ] = round( sum( co_freq.sum( 0 ) == co_freq.shape[0] ) / float( co_freq.shape[1] ), 2 )
			block_stats.loc[ block, "max_coinclusion_sub" ] = round( sum( co_freq.sum( 0 ) == co_freq.shape[0] ) / float( sum( co_freq.sum( 0 ) > 0 ) ), 2 )
		self.inclusion = pd.merge( self.inclusion, block_stats, left_index=True, right_index=True)

	def findColInclusionFreq( self ):
		self.blocks2col["block"] = self.blocks2col.index
		self.blocks2col = self.blocks2col.set_index("sample")
		self.blocks2col = self.blocks2col.loc[ self.blocks2col.block != "block", ]
		for i in self.run_composition.iterkeys():
			self.blocks2col.loc[ self.run_composition[ i ]["cols"], "r_in" ] = self.blocks2col.loc[ self.run_composition[ i ]["cols"], "r_in" ] + float( 1 )/self.nruns

	def findBlockFreq (self ):
		self.blocks["r_in"] = 0
		for i in self.run_composition.iterkeys():
			self.blocks.loc[ self.run_composition[ i ]["blocks"], "r_in" ] = self.blocks.loc[ self.run_composition[ i ]["blocks"], "r_in" ] + float( 1 )/self.nruns

	def pickCols_single( self, n ):
		"""Pick columns for a cMonkey run using predefined blocks and based on their current representation in the ensemble"""
		# what is max conditions for this run?
		n_cols = self.nruns+1
		while n_cols > self.nruns:
			n_cols = round( self.avg_col_size + ( self.avg_col_size/4 )*random.gammavariate(1,2) )
		
		# first choose excluded blocks
		excluded = self.weightedReservoirSample( 1, self.exclusion["p"] )[0]
		# update its weight
		self.exclusion.loc[ excluded, "r_out" ] = self.exclusion.loc[ excluded, "r_out" ] + 1.0/self.nruns
		self.exclusion[ "p" ] = np.power( ( 1+10e-10 ) - self.exclusion[ "r_out" ], abs( self.exclusion_percentage/100.0 - self.exclusion[ "r_out" ] ) / 1 )

		excluded_diff = set( self.blocks.index ).difference( set( excluded.split( ":::" ) ) )

		# pick a single from blocks remaining
		blocks = []
		cols = []
		while len(cols) < n_cols:
			#print blocks
			if len(cols) == 0:
				# this is the first condition
				#weights = self.blocks.loc[ excluded_diff, "p" ]
				block_one = self.weightedReservoirSample( 1, self.blocks.loc[ excluded_diff, "p" ] )[0] 
				col_one = self. translate_blocks( block_one, excluded )
				blocks.append( block_one )
				[cols.append( i ) for i in col_one if not (i in cols) ]
				[self.count_cols(i) for i in col_one if not (i in cols)]
				excluded_diff = excluded_diff.difference( set( [ block_one ] ) )
				self.updateWeights( blocks )
			else:
				#weights = self.blocks.loc[ excluded_diff, "p" ] * self.combinedWeights( excluded_diff, blocks )
				block_one = self.weightedReservoirSample( 1, self.blocks.loc[ excluded_diff, "p" ] * self.combinedWeights( excluded_diff, blocks ) )[0] 
				col_one = self. translate_blocks( block_one, excluded )
				blocks.append( block_one )
				[cols.append( i ) for i in col_one if not (i in cols) ]
				[self.count_cols(i) for i in col_one if not (i in cols)]
				excluded_diff = excluded_diff.difference( set( [ block_one ] ) )
				self.updateWeights( blocks )

		self.run_composition[n] = {
		"blocks" : blocks,
		"cols": cols,
		"excluded": excluded
		}
		return None


	def report( self, file = None ):
		
		if file == None:
			file ="./ensembleReport_"

		print "Writing reports"

		self.findBlockInclusionFreq()
		self.findColInclusionFreq()
		self.findBlockFreq()

		run_df = pd.DataFrame( 0, index = range( 1,self.nruns+1 ), columns = [ "ncols", "excluded", "blocks", "cols" ] ) 

		for i in self.run_composition.iterkeys():
			run_df.loc[ i,"ncols" ] = len( self.run_composition[ i ][ "cols" ] )
			run_df.loc[ i,"excluded" ] = self.run_composition[ i ][ "excluded" ] 
			run_df.loc[ i,"blocks" ] = self.run_composition[ i ][ "blocks" ] 
			run_df.loc[ i,"cols" ] = self.run_composition[ i ][ "cols" ] 

		run_df.to_csv( path_or_buf = open ( file + "runs.csv", mode = "w"  ), cols = [ "ncols", "excluded", "blocks", "cols" ], index_label = "run_num" )	
		self.exclusion.to_csv( path_or_buf = open ( file + "exclusionBlocks.csv", mode = "w"  ), cols = [ "block.sample.num", "r_out" ], header = [ 'block_sample_num', "excluded_freq" ], index_label = "exclusion_blocks" )
		self.blocks.to_csv( path_or_buf = open ( file + "blocks.csv", mode = "w"  ), cols = [ "block.sample.num", "r_in" ], header = [ "block_sample_num", "ensemble_freq" ], index_label = "block" )
		self.blocks2col.to_csv( path_or_buf = open ( file + "cols.csv", mode = "w"  ), cols = [ "block", "r_in" ], header = [ "block", "ensemble_freq" ], index_label = "sample" )
		self.inclusion.to_csv( path_or_buf = open ( file + "inclusionBlocks.csv", mode = "w"  ), cols = [ "block.sample.num", 'freq_single', 'freq_coinclusion', 'max_coinclusion_all', 'max_coinclusion_sub' ], header = [ "block_sample_num", 'freq_single', 'freq_coinclusion', 'max_coinclusion_all', 'max_coinclusion_sub' ], index_label = "inclusion_blocks" )

	def pickCols_all( self ):
		# get random blocks
		if self.random_blocks:
			if self.n_rand_exclusion > 0:
				self.chooseRandomBlocks()
		for i in range(1, self.nruns+1):
			if i%100 == 0:
				print str( round( float( i )/ float( self.nruns ) *100 ) ) + "% done"
			self.pickCols_single( i )
		self.report( self.report_file )
		print 'Done'

