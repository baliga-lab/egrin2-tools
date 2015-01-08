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

import pdb

import numpy as np
import pandas as pd


import sqlite3
import pymongo
from pymongo import MongoClient
import gridfs
from Bio import SeqIO
from scipy.integrate import quad

class makeCorems:

	def __init__( self, dbname = None, dbfiles = None, backbone_pval = None, out_dir = None ):

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

		if os.system( "which -s adjmat2wpairs" ) != 0:
			print "WARNING!!! You need to compile adjmat2wpairs.cpp to adjmat2wpairs and add its location to your path to detect corems"

		if os.system( "which -s compute_tanimoto" ) != 0:
			print "WARNING!!! You need to compile compute_tanimoto.cpp to compute_tanimoto and add its location to your path to detect corems"

		if os.system( "which -s cluster_communities" ) != 0:
			print "WARNING!!! You need to compile cluster_communities.cpp to cluster_communities and add its location to your path to detect corems"

		if os.system( "which -s getting_communities" ) != 0:
			print "WARNING!!! You need to compile getting_communities.cpp to getting_communities and add its location to your path to detect corems"

		if out_dir is None:
			self.out_dir = "./corem_data/"
			if not os.path.isdir( self.out_dir ):
				os.makedirs( self.out_dir )
		else:
			self.out_dir = self.out_dir
			if not os.path.isdir( self.out_dir ):
				os.makedirs( self.out_dir )
		print "Corem data will be output to:", self.out_dir

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

		def structureRowRow ( key_row, sub_row, data_counts, data_counts_norm, backbone_pval, row_row_collection ):
			try:
				weight = row_row_collection.find_one( { "row_ids" : [ self.row2id[ sub_row ], self.row2id[ key_row ] ]  } )["weight"]
				# check to see if this pair already exists and if current weight is greater
				# if current weight is greater and backbone_pval is significant, update the weight in MongoDB
				if ( data_counts_norm > weight ) & ( backbone_pval <= self.backbone_pval ):
					row_row_collection.update({ "row_ids" : [ self.row2id[ sub_row ], self.row2id[ key_row ] ] }, { "$set":{ "weight" : data_counts_norm, "backbone_pval" : backbone_pval } } )
				d = None
			except Exception:
				d = {
				"row_ids": [ self.row2id[ key_row ], self.row2id[ sub_row ] ], 
				"counts": data_counts,
				"weight": data_counts_norm,
				"backbone_pval": backbone_pval
				}
			return d

		row_row_collection = self.db.row_row
		counter = 1
		
		for i in self.row2id.keys():
			print i
			if counter%250 == 0:
				print "%s percent done" % str( round( float( counter ) / len( self.row2id.keys() ), 2 )*100 )
			# check if already exists in DB
			data_counts = self.getRowCo( i )
			# set self counts to 0 and normalize other counts
			data_counts[ i ] = 0
			data_counts = data_counts[ data_counts>0 ]
			data_counts_norm = data_counts.copy()
			data_counts_norm = data_counts_norm/ sum( data_counts_norm )
			# only keep values > 0
			backbone_data_counts = self.extractBackbone( data_counts_norm )
			
			to_write = [ structureRowRow( i, j, data_counts[j], data_counts_norm[j], backbone_data_counts[j], row_row_collection ) for j in data_counts.index ]
			to_write = filter( None, to_write)
			row_row_collection.insert( to_write )

			counter = counter + 1

	def writeEdgeList( self ):
		return None

	def runCoremCscripts( self ):
		return None




