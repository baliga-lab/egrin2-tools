#!/usr/bin/env python

"""Do all steps to assemble a cMonkey2 ensemble.

Example:

python assembler.py --organism eco --ratios /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz --targetdir /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ --ncbi_code 511145 --ensembledir /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/ --col_annot /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0

python assembler.py --organism eco --ratios ./ratios_eco_m3d.tsv.gz --targetdir ./ --ncbi_code 511145 --ensembledir ./eco-ens-m3d/ --col_annot ./E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0

python assembler.py --organism mtu --ratios ./20141130.MTB.all.ratios.csv.gz --targetdir ./ --ncbi_code 83332 --ensembledir ./  --n_resamples 0


"""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import argparse
import os
import itertools

from assemble.sql2mongoDB import *
from assemble.makeCorems import * 
from assemble.resample import *

DESCRIPTION = """assemble.py - prepare cluster runs"""

if __name__ == '__main__':

	import argparse
	import os
	import itertools

	parser = argparse.ArgumentParser(description=DESCRIPTION)
	parser.add_argument('--organism', required=True, type=str, help="Organism code")
	parser.add_argument('--ratios', required=True, help="These should be original 'raw' normalized ratios, not the standardized ratios used in cMonkey")
	parser.add_argument('--targetdir', required=True, help="Where should MongoDB and corem data be stored")
	parser.add_argument('--ncbi_code', required=True)
	parser.add_argument('--cores', default=3, type=int)
	parser.add_argument('--ensembledir', default=None, help="Location of the ensemble runs. Default: cwd")
	parser.add_argument('--col_annot', default=None, help="A tab-delimited file with condition annotations")
	parser.add_argument('--host', default="localhost", help="Host for MongoDB")
	parser.add_argument('--port', default=27017, help="MongoDB port", type=int)
	parser.add_argument('--prefix', default=None, help="Ensemble run name prefix. Default: *organism*-out-")
	parser.add_argument('--row_annot', default=None, help="Optional row (gene) annotation tab-delimited file. If not specified, annotations will be downloaded from MicrobesOnline using --ncbi_code.")
	parser.add_argument('--row_annot_matchCol', default=None, help="Name of column in row_annot that matches row names in ratios file.")
	parser.add_argument('--gre2motif', default=None, help="Motif->GRE clustering file")
	parser.add_argument('--db', default=None, help="Optional ensemble MongoDB database name")
	parser.add_argument('--genome_annot', default=None, help="Optional genome annotation file. Automatically downloaded from MicrobesOnline using --ncbi_code")
	parser.add_argument('--backbone_pval', default = 0.05, type=float )
	parser.add_argument('--link_comm_score', default = None )
	parser.add_argument('--link_comm_increment', default = None )
	parser.add_argument('--link_comm_density_score', default = None )
	parser.add_argument('--corem_size_threshold', default = None )
	parser.add_argument('--n_resamples', default=10000, type=int, help="# resamples to compute for corem condition assignment")

	args = parser.parse_args()

	# Initialize to find problems early!!
	sql2mongo = sql2mongoDB( organism = args.organism, host = args.host, port = args.port, ensembledir = args.ensembledir, prefix = args.prefix, ratios_raw = args.ratios, gre2motif = args.gre2motif, col_annot = args.col_annot, ncbi_code = args.ncbi_code, dbname = args.db, db_run_override = None, genome_file = args.genome_annot, row_annot = args.row_annot, row_annot_match_col = args.row_annot_matchCol )
	if len( sql2mongo.db_files ) >0:
		# Merge sql into mongoDB
		sql2mongo.compile()
		corems = makeCorems( organism = args.organism, host = args.host, port = args.port, db = args.db, dbfiles = None, backbone_pval = args.backbone_pval, out_dir = args.targetdir, n_subs = args.cores, link_comm_score = args.link_comm_score, link_comm_increment = args.link_comm_increment, link_comm_density_score = args.link_comm_density_score, corem_size_threshold = args.corem_size_threshold )

		# Make corems
		corems.rowRow()
		corems.runCoremCscripts()
		corems.getCorems()
		corems.addCorems()

		if args.n_resamples > 0:
			# Make resample database
			print "Computing resamples"
			client = MongoClient( host = args.host, port= args.port )
			db = sql2mongo.dbname 
			cols = range( 0,client[ db ][ "col_info" ].count( ) )
			corem_sizes = list( set( [ len( i[ "rows" ] ) for i in client[ db ][ "corem" ].find( {}, {"rows":1} ) ] ) )
			corem_sizes.sort( )
			tmp = Parallel(n_jobs=args.cores )( delayed( colResampleInd )( args.host, sql2mongo.dbname, i, cols, n_resamples = args.n_resamples) for i in corem_sizes )

			# Finish corems by adding computed resamples
			print "Finishing corems"
			print "Adding condition information"
			corems.finishCorems()

		outfile =  sql2mongo.prefix + str(datetime.datetime.utcnow()).split(" ")[0] + ".mongodump"
		print "Writing EGRIN2 MongoDB to %s" % sql2mongo.targetdir + outfile  
		sql2mongo.mongoDump( sql2mongo.dbname, outfile )

		print "Done"
	else:	
		print "No cMonkey2 DB files. Please specify --ensembledir"
		