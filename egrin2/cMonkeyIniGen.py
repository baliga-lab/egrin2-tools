#!/usr/bin/env python

"""
cMonkey ini generator for customizing ensemble runs
"""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import os.path
import random

class cMonkeyIniGen:     

	def __init__( self, params = {} ):
		# params is a dictionary of parameters 

		# [General]
		if "normalize_ratios" not in params.keys():
			params[ "normalize_ratios" ] = True
		if "num_iterations" not in params.keys():
			params[ "num_iterations" ] = 2000
		if "start_iteration" not in params.keys():
			params[ "start_iteration" ] = 1
		if "output_dir" not in params.keys():
			params[ "output_dir" ] = "out"
		if "cache_dir" not in params.keys():
			params[ "cache_dir" ] = "cache"
		if "tmp_dir" not in params.keys():
			params[ "tmp_dir" ] = ""
		if "dbfile_name" not in params.keys():
			params[ "dbfile_name" ] = "cmonkey_run.db"
		if "use_mutliprocessing" not in params.keys():
			params[ "use_multiprocessing" ] = True
		if "checkpoint_interval" not in params.keys():
			params[ "checkpoint_interval" ] = 99999
		if "num_core" not in params.keys():
			params[ "num_cores" ] = 1
		if "stats_frequency" not in params.keys():
			params[ "stats_frequency" ] = 10
		if "results_frequency" not in params.keys():
			params[ "result_frequency"] = 10
		if "debug_frequency" not in params.keys():
			params[ "debug_frequency"] = 50
		if "postadjust" not in params.keys():
			params[ "postadjust" ] = True
		if "add_fuzz" not in params.keys():
			params[ "add_fuzz" ] = "rows"
		if "num_clusters" not in params.keys():
			params[ "num_clusters" ] = random.randint(150, 550)
		if "random_seed" not in params.keys():
			params[ "random_seed" ] = 1
			random.seed( random_seed )
		if "log_subresults" not in params.keys():
			params[ "log_subresults" ] = True
		if "case_sensitive" not in params.keys():
			params[ "case_sensitive" ] = True
		if "rsat_base_url" not in params.keys():
			params[ "rsat_base_url" ] = "http://embnet.ccg.unam.mx/rsa-tools"
		if "organism_code" not in params.keys():
			params[ "organism_code" ] = ""
		if "use_operons" not in params.keys():
			params[ "use_operons" ] = random.randint(0, 1) == 1
		if "use_string" not in params.keys():
			params[ "use_string" ] = random.randint(0, 1) == 1
		if "checkratios" not in params.keys():
			params[ "checkratios" ] = False

		# [Membership]
		if "probability_row_change" not in params.keys():
			params[ "probability_row_change" ] = 0.5
		if "probability_column_change" not in params.keys():
			params[ "probability_column_change" ] = 1.0
		if "max_changes_per_row" not in params.keys():
			params[ "max_changes_per_row" ] = 1
		if "max_changes_per_column" not in params.keys():
			params[ "max_changes_per_column" ] = 5
		if "min_cluster_rows_allowed" not in params.keys():
			params[ "min_cluster_rows_allowed" ] = 3
		if "max_cluster_rows_allowed" not in params.keys():
			params[ "max_cluster_rows_allowed" ] = 70
		if "clusters_per_row" not in params.keys():
			params[ "clusters_per_row" ] = 2
		if "clusters_per_column" not in params.keys():
			params[ "clusters_per_column" ] = 215

		# [Scoring]
		if "quantile_normalize" not in params.keys():
			params[ "quantile_normalize" ] = False

		# [Rows]
		if "row_schedule" not in params.keys():
			params[ "row_schedule" ] = "1,2"
		if "row_scaling_const" not in params.keys():
			params[ "row_scaling_const" ] = random.uniform(0.0, 3.0)

		# [Columns]
		if "col_schedule" not in params.keys():
			params[ "col_schedule" ] = "1,5"

		# [Networks]
		if "networks_schedule" not in params.keys():
			params[ "networks_schedule" ] = "1,7"
		if "networks_scaling_rvec" not in params.keys():
			params[ "networks_scaling_rvec" ]  = "seq(1e-5, 0.5, length=num_iterations*3/4)"
		if "networks_scaling_const" not in params.keys():
			params[ "networks_scaling_const" ] = random.uniform(0.0, 1.0)

		# [Motifs]
		if "sequence_types" not in params.keys():
			params[ "sequence_types" ] = "upstream"
		if "motifs_schedule" not in params.keys():
			params[ "motifs_schedule" ] = "2,10"
		if "motifs_scaling_rvec" not in params.keys():
			params[ "motifs_scaling_rvec" ] = "c(rep(1e-5, 100), seq(1e-5, 1, length=num_iterations*3/4))"

		# [MEME]
		BGORDER = [None, 0, 1, 2, 3, 4, 5]

		if "meme_global_background" not in params.keys():
			params[ "meme_global_background" ] = True
		if "meme_schedule" not in params.keys():
			params[ "meme_schedule" ] = "1,100"
		if "meme_nmotifs_rvec" not in params.keys():
			params[ "meme_nmotifs_rvec" ] = "c(rep(1, num_iterations/3), rep(2, num_iterations/3))"
		if "use_revcomp" not in params.keys():
			params[ "use_revcomp" ] = True
		if "max_width" not in params.keys():
			params[ "max_width" ] = random.randint(12, 30)
		if "background_order" not in params.keys():
			params[ "background_order" ] = BGORDER[ random.randint( 0, 6 ) ]
		if "arg_mod" not in params.keys():
			params[ "arg_mod" ] = "zoops"
		#
		# I don't know what the hell these do, but they were in Wei-ju's code!!!
		#
		if "meme_string_weight" not in params.keys():
			params[ "meme_string_weight" ] = random.uniform(0.0, 1.0) * 0.5 + 0.2
        		if "meme_operon_weight" not in params.keys():
        			params[ "meme_operon_weight" ] = random.uniform(0.0, 1.0) * 0.5 + 0.2

		# [Weeder]
		if "weeder_global_background" not in params.keys():
			params[ "weeder_global_background" ] = True
		if "weeder_schedule" not in params.keys():
			params[ "weeder_schedule" ] = "1,100"
		if "weeder_nmotifs_rvec" not in params.keys():
			params[ "weeder_nmotifs_rvec" ] = "c(rep(1, num_iterations/3), rep(2, num_iterations/3))"
		if "orgcode" not in params.keys():
			params[ "orgcode" ] = ""
		if "freqfile_dir" not in params.keys():
			params[ "freqfile_dir" ] = ""
		if "analysis" not in params.keys():
			params[ "analysis" ] = "small"
		if "top" not in params.keys():
			params[ "top" ] = 50

		# [SequenceType-upstream]
		if "search_distance" not in params.keys():
			params[ "search_distance" ] = str( random.randint( -20, 0 ) ) + "," + str( random.randint( 100, 200 ) )
		if "scan_distance" not in params.keys():
			params[ "scan_distance" ] = str( random.randint( -50, 0 ) ) + "," + str( random.randint( 150, 250 ) )

		# [SetEnrichment]
		if "setenrichment_schedule" not in params.keys():
			params[ "setenrichment_schedule" ] = "1,7"
		if "setenrichment_scaling_rvec" not in params.keys():
			params[ "setenrichment_scaling_rvec" ] = "seq(1e-5, 1.5, length=num_iterations*3/4)"
		if "set_types" not in params.keys():
			# can be array!
			params[ "set_types" ] = ""

		# [SetEnrichment-set_types]
		if "set_file" not in params.keys():
			# can be an array!
			params[ "set_file" ] = ""
		if "set_weight" not in params.keys():
			# can be an array!
			params[ "set_weight" ] = ""

		self.CMONKEY_INI_TEMPLATE =  """
		[General]
		normalize_ratios = %(normalize_ratios)s
		num_iterations = %(num_iterations)s
		start_iteration = %(start_iteration)s
		output_dir = %(output_dir)s
		cache_dir = %(cache_dir)s
		tmp_dir = %(tmp_dir)s
		dbfile_name = %(dbfile_name)s
		use_multiprocessing = %(use_multiprocessing)s
		checkpoint_interval = %(checkpoint_interval)s
		num_cores = %(num_cores)s
		stats_frequency = %(stats_frequency)s
		result_frequency = %(result_frequency)s
		debug_frequency = %(debug_frequency)s
		postadjust = %(postadjust)s
		add_fuzz = %(add_fuzz)s
		num_clusters = %(num_clusters)s
		random_seed = %(random_seed)s
		log_subresults = %(log_subresults)s
		case_sensitive = %(case_sensitive)s
		rsat_base_url = %(rsat_base_url)s
		organism_code = %(organism_code)s
		use_operons = %(use_operons)s
		use_string = %(use_string)s
		checkratios = %(checkratios)s

		[Membership]
		probability_row_change = %(probability_row_change)s
		probability_column_change = %(probability_column_change)s
		max_changes_per_row = %(max_changes_per_row)s
		max_changes_per_column = %(max_changes_per_column)s
		min_cluster_rows_allowed = %(min_cluster_rows_allowed)s
		max_cluster_rows_allowed = %(max_cluster_rows_allowed)s
		clusters_per_row = %(clusters_per_row)s
		clusters_per_column = %(clusters_per_column)s

		[Scoring]
		quantile_normalize = %(quantile_normalize)s

		[Rows]
		schedule = %(row_schedule)s
		scaling_const = %(row_scaling_const)s

		[Columns]
		schedule = %(col_schedule)s

		[Networks]
		schedule = %(networks_schedule)s
		scaling_rvec = %(networks_scaling_rvec)s
		scaling_const = %(networks_scaling_const)s

		[Motifs]
		sequence_types = %(sequence_types)s
		schedule = %(motifs_schedule)s
		scaling_rvec = %(motifs_scaling_rvec)s

		[MEME]
		global_background = %(meme_global_background)s
		schedule = %(meme_schedule)s
		nmotifs_rvec = %(meme_nmotifs_rvec)s
		use_revcomp = %(use_revcomp)s
		max_width = %(max_width)s
		background_order = %(background_order)s
		arg_mod = %(arg_mod)s
		string_weight = %(meme_string_weight)s
		operon_weight = %(meme_operon_weight)s

		[Weeder]
		global_background = %(weeder_global_background)s
		schedule = %(weeder_schedule)s
		nmotifs_rvec = %(weeder_nmotifs_rvec)s
		orgcode = %(orgcode)s
		freqfile_dir = %(freqfile_dir)s
		analysis = %(analysis)s
		top = %(top)s

		[SequenceType-upstream]
		search_distance = %(search_distance)s
		scan_distance = %(scan_distance)s

		[SetEnrichment]
		schedule = %(setenrichment_schedule)s
		scaling_rvec= %(setenrichment_scaling_rvec)s
		set_types = %(set_types)s""" % params

		for i in range(0, len( params["set_types"].split( "," ) ) ):
			if params["set_types"] != "":
				set_types = params["set_types"].split( "," )[i]
				set_file = params[ "set_file" ].split( "," )[i]
				set_weight = params[ "set_weight" ].split( "," )[i]
				SETENRICHMENT_TEMPLATE = """
				[SetEnrichment-%(set_types)s]
				set_file = %(set_file)s
				weight = %(set_weight)s
				"""
				self.CMONKEY_INI_TEMPLATE = self.CMONKEY_INI_TEMPLATE + SETENRICHMENT_TEMPLATE

