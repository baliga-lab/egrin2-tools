#!/usr/bin/env python

"""Generate corems using cMoneky ensemble MongoDB. Add them to an existing MongoDB"""

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

	def __init__( self, dbname = None, dbfiles = None, backbone_pval = None ):

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

		if backbone_pval == None:
			self.backbone_pval = 0.05
		else:
			self.backbone_pval = backbone_pval

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
			if pval <= 0.05:
				backbone_data_counts[i] = pval
			else:
				backbone_data_counts[i] = 0
		return backbone_data_counts


	def rowRow( self ):
		"""Construct row-row co-occurrence matrix (ie gene-gene co-occurrence)"""

		def structureRowRow ( key_row, sub_row, data_counts, data_counts_norm, backbone_data_counts ):
			d = {
			"row_ids": [ self.row2id[ i ], self.row2id[ j ] ], 
			""
			}


		row_row_collection = self.db.row_row
		counter = 1
		
		for i in self.row2id.keys():
			if counter%250 == 0:
				print "%s percent done" % str( round( float( counter ) / len( self.row2id.keys )*100, 1 ) )
			# check if already exists in DB
			data_counts = self.getRowCo( i )
			# set self counts to 0 and normalize other counts
			data_counts[ i ] = 0
			data_counts = data_counts[ data_counts>0 ]
			data_counts_norm = data_counts.copy()
			data_counts_norm = data_counts_norm/ sum( data_counts_norm )
			# only keep values > 0
			backbone_data_counts = self.extractBackbone( data_counts_norm )
			
			to_write = [ structureRowRow( i, j, data_counts, data_counts_norm, backbone_data_counts ) for j in data_counts.index ]
			col_info_collection.insert( d_f )

			counter = counter + 1

	def runCoremCscripts( self ):





