#!/usr/bin/env python

"""Initialize  mongo database from individual cMonkey 
runs (SQLite) plus some additional tables and run_id column, 
including motifs and motif clusters"""

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

import numpy as np

import sqlite3
from pymongo import MongoClient

# connect to database
# make sure mondod is running
retvalue = os.system("nohup mongod --port 27017 --quiet &")

client = MongoClient('mongodb://localhost:27017/') 

db = client.test_database

collection = db.test_collection

class sql2mongoDB:
    
    def __init__( self, e_dir = None, prefix = None, gene_info = None ):
    	if prefix is None:
    		prefix = 'eco-out-'
    	else:
        		self.prefix = prefix
        	if e_dir is None:
    		e_dir = '/eco-ens-m3d/'
    	else:
        		self.prefix = prefix
        	self.files = np.sort( np.array( glob.glob( e_dir + prefix + "???/cmonkey_run.db" ) ) ) # get all cmonkey_run.db files

    def parseSQLite( self, db_file ):
    	"""Create python dictionaries for bulk import into MongoDB collections"""
        	conn = sqlite3.connect(db_file)
        	c = conn.cursor()

