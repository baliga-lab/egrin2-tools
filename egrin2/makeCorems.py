#!/usr/bin/env python

"""Generate corems using cMonkey ensemble MongoDB. Add them to an existing MongoDB"""

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
from zipfile import ZipFile
from itertools import combinations
import subprocess
import shutil
import decimal
from math import ceil

import pdb

import numpy as np
import pandas as pd


import sqlite3
import pymongo
from pymongo import MongoClient
from Bio import SeqIO
from scipy.integrate import quad
import matplotlib.pyplot as plt

class makeCorems:

	def __init__( self, dbname = None, dbfiles = None, backbone_pval = None, out_dir = None, n_subs = None, link_comm_score = None, link_comm_increment = None, link_common_density_score = None, corem_size_threshold = None ):

		# connect to database
		# make sure mongodb is running
		retvalue = os.system("nohup mongod --port 27017 --quiet &")
		
		try:
			client = MongoClient('mongodb://localhost:27017/') 
			print "Connected to MongoDB"
		except pymongo.errors.ConnectionFailure, e:
			print "Could not connect to MongoDB: %s" % e
			return None

		if dbname == None:
			self.dbname = "egrin2_db"
		else:
			self.dbname = dbname

		if self.dbname in client.database_names():
			print "Found ensemble database: %s" % self.dbname
		else:
			print "Could not locate a MongoDB database with name: %s\nAttempting to perform DB import at: %s" % ( self.dbname, dbfiles )
			if dbfiles != None:
				try:
					self.mongoRestore( self.dbname, dbfiles )
				except Exception:
					return None
			else:
				return None
		
		self.db = client[self.dbname]
		
		self.row2id = {}
		self.id2row = {}
		for i in self.db.row_info.find( { }, { "egrin2_row_name": "1", "row_id": "1" } ):
			self.row2id[ i[ "egrin2_row_name" ] ] = i[ "row_id" ]
			self.id2row[ i[ "row_id" ] ] = i[ "egrin2_row_name" ]

		self.col2id =  {}
		self.id2col =  {}
		for i in self.db.col_info.find( { }, { "egrin2_col_name": "1", "col_id": "1" } ):
			self.col2id[ i[ "egrin2_col_name" ] ] = i[ "col_id" ]
			self.id2col[ i[ "col_id" ] ] = i[ "egrin2_col_name" ]

		if backbone_pval is None:
			self.backbone_pval = 0.05
		else:
			self.backbone_pval = backbone_pval

		self.cFail = False
		if os.system( "which -s adjmat2wpairs" ) != 0:
			print "WARNING!!! You need to compile adjmat2wpairs.cpp to adjmat2wpairs and add its location to your path to detect corems"
			self.cFail = True

		if os.system( "which -s compute_tanimoto" ) != 0:
			print "WARNING!!! You need to compile compute_tanimoto.cpp to compute_tanimoto and add its location to your path to detect corems"
			self.cFail = True

		if os.system( "which -s cluster_communities" ) != 0:
			print "WARNING!!! You need to compile cluster_communities.cpp to cluster_communities and add its location to your path to detect corems"
			self.cFail = True

		if os.system( "which -s getting_communities" ) != 0:
			print "WARNING!!! You need to compile getting_communities.cpp to getting_communities and add its location to your path to detect corems"
			self.cFail = True

		if out_dir is None:
			self.out_dir = "corem_data"
			if not os.path.isdir( self.out_dir ):
				os.makedirs( self.out_dir )
		else:
			self.out_dir = self.out_dir
			if not os.path.isdir( self.out_dir ):
				os.makedirs( self.out_dir )
		print "Corem data will be output to:", self.out_dir

		if n_subs is None:
			# number of subprocesses to spawn
			self.n_subs = 4

		if link_comm_score is None:
			# use link similarity def of (0) Ahn or (1) Kalinka
			self.link_comm_score = 0
		else:
			self.link_comm_score = link_comm_score 

		if link_comm_increment is None: 
			# amount to increment community detection threshold
			self.link_comm_increment = 0.1
		else:
			self.link_comm_increment = link_comm_increment

		if link_common_density_score is None:
			# score used to evaluate global density of communities (1,2,3,4,5) 
			self.link_common_density_score = 5
		else:
			self.link_common_density_score = link_common_density_score

		if corem_size_threshold is None:
			# minimum size of corem, # edges
			self.corem_size_threshold = 3
		else: 
			self.corem_size_threshold = corem_size_threshold

	def mongoRestore( self, db, infile ):
		"""Read contents of binary MongoDB dump into MongoDB instance"""
		sys_command = "mongorestore --db " + db + " " + infile
		print sys_command
		os.system(sys_command)

	def getRowCo( self, row_id ):
		"""Given a row (gene), count all of the other rows that occur with it in a bicluster"""
		data = []
		for i in self.db.bicluster_info.find( { "rows" : { "$all" : [ self.row2id[ row_id ] ] } }, { "rows": "1" } ):
			for j in i[ "rows" ]:
				try:
					data.append( self.id2row[ j ] )
				except:
					continue
		data_counts = pd.Series( data ).value_counts()
		return data_counts

	def extractBackbone( self, data_counts ):
		"""Extract the significant elements from rBr co-occurrence matrix"""
		backbone_data_counts = data_counts.copy()

		def integrand( x, k ):
			 return np.power( 1.0-x , k-2.0 )

		for i in data_counts.index:
			k = len( data_counts )
			pval = 1-(k-1)*quad( integrand, 0, data_counts[ i ], args=( k ) ) [0]
			backbone_data_counts[i] = pval
		return backbone_data_counts

	def rowRow( self ):
		"""Construct row-row co-occurrence matrix (ie gene-gene co-occurrence)"""
		
		row_row_collection = self.db.row_row

		# remove existing edgeList file if it exists
		if os.path.exists( os.path.abspath( os.path.join( self.out_dir,"edgeList" ) ) ):
			print "Found edgeList file at %s. Removing it." % os.path.abspath( os.path.join( self.out_dir,"edgeList" ) )
			os.remove ( os.path.abspath( os.path.join( self.out_dir,"edgeList" ) ) )
			print "Dropping row_row MongoDB collection as a precaution."
			row_row_collection.drop()

		def addToD( d, ind1, ind2, val ):
			if ind1 not in d.iterkeys():
				d[ ind1 ] = {}
			d[ ind1 ][ ind2 ] = val
			return d

		def writeRowRow( row ):
			"""Only write rows with significant backbone pvals to edgeList file"""
			if float( row["backbone_pval"] ) <= self.backbone_pval:
				f.write( ( " " ).join( [ self.id2row[ row[ "row_ids" ][ 0 ] ], self.id2row[ row[ "row_ids" ][ 1 ] ], str( row[ "weight" ] ), "\n" ] ) )

		def structureRowRow ( key_row, sub_row, data_counts, data_counts_norm, backbone_pval, row_row_collection ):
			try:
				# check to see if this pair already exists and if current weight is greater
				weight = self.rowrow_ref[ self.row2id[ sub_row ] ][ self.row2id[ key_row ] ]
				
				if ( data_counts_norm > weight ) & ( backbone_pval <= self.backbone_pval ):
					# if current weight is greater and backbone_pval is significant, update the weight in MongoDB
					row_row_collection.update({ "row_ids" : [ self.row2id[ sub_row ], self.row2id[ key_row ] ] }, { "$set":{ "weight" : data_counts_norm, "backbone_pval" : backbone_pval } } )
					self.rowrow_ref = addToD( self.rowrow_ref, self.row2id[ sub_row ], self.row2id[ key_row ], data_counts_norm )
				d = None
			except Exception:
				self.rowrow_ref = addToD( self.rowrow_ref, self.row2id[ sub_row ], self.row2id[ key_row ], data_counts_norm )
				d = {
				"row_ids": [ self.row2id[ key_row ], self.row2id[ sub_row ] ], 
				"counts": data_counts,
				"weight": data_counts_norm,
				"backbone_pval": backbone_pval
				}
			return d

		counter = 1
		
		# make a dictionary to keep track of rows in db and their weights
		# was too slow using mongoDB lookups...
		self.rowrow_ref = {}
		print "Constructing row-row co-occurrence matrix. This will take some time..."
		for i in self.row2id.keys():
			#print i
			if counter%250 == 0:
				print "%s percent done" % str( round( float( counter ) / len( self.row2id.keys() ), 2 )*100 )
			# check if already exists in DB
			data_counts = self.getRowCo( i )
			# set self counts to 0 and normalize other counts
			data_counts[ i ] = 0
			data_counts = data_counts[ data_counts>0 ]
			data_counts_norm = data_counts.copy()
			data_counts_norm = data_counts_norm / sum( data_counts_norm )
			# only keep values > 0
			backbone_data_counts = self.extractBackbone( data_counts_norm )
			
			to_write = [ structureRowRow( i, j, data_counts[j], data_counts_norm[j], backbone_data_counts[j], row_row_collection ) for j in data_counts.index ]
			to_write = filter( None, to_write)

			# write edgeList file
			with open( os.path.abspath( os.path.join( self.out_dir,"edgeList" ) ), mode = "a+") as f:
				  [  writeRowRow( j ) for j in to_write ]

			row_row_collection.insert( to_write )

			counter = counter + 1
		
		# clean up
		del self.rowrow_ref

	def runCoremCscripts( self ):

		def drange(start, stop = None, step = 1, precision = None):
			"""drange generates a set of Decimal values over the
			range [start, stop) with step size step

			drange([start,] stop, [step [,precision]])"""

			try:
			    _xrange = xrange
			except NameError:
			    _xrange = range

			if stop is None:
			    for x in _xrange(int(ceil(start))):
			        yield x
			else:
			    # find precision
			    if precision is not None:
			        decimal.getcontext().prec = precision
			    # convert values to decimals
			    start = decimal.Decimal(start)
			    stop = decimal.Decimal(stop)
			    step = decimal.Decimal(step)
			    # create a generator expression for the index values
			    indices = (
			        i for i in _xrange(
			            0, 
			            ((stop-start)/step).to_integral_value()
			        )
			    )  
			    # yield results
			    for i in indices:
			        yield float(start + step*i)
	
		if self.cFail:
			print "Cannot detect corems because one or more community detection C++ scrips are either (1) not compiled, (2) not in the $PATH variable, or (3) incorrectly named. Resolve previous warning."
			return None
		else:
			p = subprocess.Popen( [ "adjmat2wpairs", "edgeList", "0", "0" ], cwd=os.path.abspath( self.out_dir ) )
			p.wait()
			
			ranges = range( 0, len(self.row2id)+1, ( len( self.row2id ) + 1 ) / self.n_subs )
			# make sure last is # genes
			ranges[len(ranges)-1] = len(self.row2id) + 1
			# "0" parameter is for link similarity of Ahn et al. Use "1" for measure proposed by Kalinka et al
			command_template = [ 'compute_tanimoto', 'edgeList', str( self.link_comm_score ), '0', len( self.row2id ) ]
			commands = []
			for i in range( 0, len( ranges ) - 1):
				cmd = list( command_template )
				cmd[ 3 ] = str( ranges[ i ] )
				cmd[ 4 ] = str( ranges[ i+1] )
				commands.append( cmd  ) 
			# run in parallel
			processes = [ subprocess.Popen( cmd, cwd=os.path.abspath( self.out_dir ) ) for cmd in commands ]

			# wait for completion
			for p in processes: p.wait()

			# clean up
			with open( os.path.join( os.path.abspath( self.out_dir ),"edgeList.tanimoto" ), 'w' ) as outfile:
			    for infile in glob.glob( os.path.join( os.path.abspath( self.out_dir ),"edgeList.tanimoto_*" ) ):
			        shutil.copyfileobj( open( infile ), outfile )
			p = subprocess.Popen( [ "rm edgeList.tanimoto_*" ], cwd= os.path.abspath( self.out_dir ), shell = True )
			p.wait()

			print "Clustering link communities across thresholds defined by increment:", self.link_comm_increment
			command_template = [ 'cluster_communities', 'edgeList', 1]
			commands = []
			for i in list( drange( self.link_comm_increment,1 + self.link_comm_increment, self.link_comm_increment, 1 ) ):
				cmd = list( command_template )
				cmd[ 2 ] = str( i )
				commands.append( cmd  ) 
			# run in parallel
			count_1 = 0
			while count_1 < len( commands ):
				# only use define number of subprocesses at a time
				count_2 = count_1 + self.n_subs
				processes = [ subprocess.Popen( cmd, cwd=os.path.abspath( self.out_dir ) ) for cmd in commands[ count_1:count_2 ] ]
				# wait for completion
				for p in processes: p.wait()
				count_1 = count_2 
				if count_2 < len( commands ):
					print "%s percent done" % ( round( float( count_2 ) / len( commands ), 2)*100 )
				
			# clean up
			with open( os.path.join( os.path.abspath( self.out_dir ),"edgeList.density" ), 'w' ) as outfile:
			    for infile in glob.glob( os.path.join( os.path.abspath( self.out_dir ),"edgeList.density_*" ) ):
			        shutil.copyfileobj( open( infile ), outfile )
			p = subprocess.Popen( [ "rm edgeList.density_*" ], cwd= os.path.abspath( self.out_dir ), shell = True )
			p.wait()

			# choose cutoff
			density = pd.read_csv( os.path.join( os.path.abspath( self.out_dir ),"edgeList.density" ), sep = "\t", header = None )
			# map density score to column
			score_map = { 1:2, 2:4, 5:6, 3:7, 4:8 }
			# plot

			plt.plot(density.loc[:,0], density.loc[:,2], 'r-')
			
			plt.axis([0, 1, 0, 2])
			plt.show()

			fig, ax1 = plt.subplots()
			t = np.arange(0.01, 10.0, 0.01)
			s1 = np.exp(t)
			ax1.plot(t, s1, 'b-')
			ax1.set_xlabel('time (s)')
			# Make the y-axis label and tick labels match the line color.
			ax1.set_ylabel('exp', color='b')
			for tl in ax1.get_yticklabels():
			    tl.set_color('b')


			ax2 = ax1.twinx()
			s2 = np.sin(2*np.pi*t)
			ax2.plot(t, s2, 'r.')
			ax2.set_ylabel('sin', color='r')
			for tl in ax2.get_yticklabels():
			    tl.set_color('r')
			plt.show()


			
			






