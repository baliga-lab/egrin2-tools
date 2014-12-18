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

	def __init__( self, dbname = None, dbfiles = None ):

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

	def extractBackbone( self, row ):
		"""Extract the significant elements from rBr co-occurrence matrix"""
		def integrand( x, k ):
			return ( 1-x )^( k-2 )
  		# fcn to calculate degree,k
  		def calc_k( i ):
  			sum( i > 0 )
		row_v = self.rBr.loc[ row, ]
		for i in tmp:
			if i > 0:
				1-(k-1)*sci.integrate.quad( integrand, 0, i, args=( k ) ) 
			else:
				return 0

	def rowRow( self ):
		"""Construct gene-gene co-occurence matrix"""
		self.rBr = pd.DataFrame( 0, index = self.row2id.keys(), columns = self.row2id.keys() )
		counter = 1
		for i in self.rBr.index:
			if counter%250 == 0:
				print "%s percent done" % str( round( float( counter ) / len( self.rBr.index )*100, 1 ) )
			data_counts = self.getRowCo( i  )
			self.rBr.loc[ i, data_counts.index ] = data_counts.values
			counter = counter + 1




