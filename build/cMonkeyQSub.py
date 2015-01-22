#!/usr/bin/env python
"""cMonkeyQSub.py
This is a tool to prepare a directory for ensemble runs. Given an organism code
and the ratios, matrix, split up the matrix into sub matrices and create
SGE qsub shell scripts

Example:

With pre-defined condition blocks:

python cMonkeyQSub.py --organism mtu --ratios 20141130.MTB.all.ratios.csv --targetdir mtu-ens-20140120 --numruns 500 --blocks 20141202.MTB.EGRIN2.blocks.csv --inclusion 20141202.MTB.EGRIN2.inclusion.blocks.csv --exclusion 20141202.MTB.EGRIN2.exclusion.blocks.csv --pipeline setenrich_pipeline.json --setenrich chipseq,tfoe --setenrich_files 20140725.MTB.ChIPSeq.csv,20140725.MTB.DE.csv --csh

Without pre-defined condition blocks:

python cMonkeyQSub.py --organism mtu --ratios 20141130.MTB.all.ratios.csv --targetdir mtb-ens-20141230 --numruns 20 --mincols 50 --num_cores 1 --csh
"""

__author__ = "Aaron Brooks, Wei-ju Wu"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks", "Wei-ju Wu"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import argparse
import os
import itertools

#import cmonkey.datamatrix as dm
# need to be in python path!!!
from cMonkeyIniGen import *
from ensemblePicker import * 

DESCRIPTION = """ensemble.py - prepare cluster runs"""

# Templates for Bourne Shell
QSUB_TEMPLATE_HEADER = """#!/bin/bash

export LD_LIBRARY_PATH=/tools/lib:/tools/R-3.0.3/lib64/R/lib
export PATH=/tools/bin:${PATH}
export BATCHNUM=`printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE = """#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -t 1-%d
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

python cmonkey.py --organism %s --ratios %s --config %s --out %s

bzip2 -f %s/*.pkl
"""

# Templates for csh

QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh

setenv LD_LIBRARY_PATH /tools/lib:/tools/R-3.0.3/lib64/R/lib
setenv PATH /tools/bin:${PATH}
setenv BATCHNUM `printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE_CSH = """#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -t 1-%d
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

python cmonkey.py --organism %s --ratios %s --config %s --out %s --minmize_io

bzip2 -f %s/*.pkl
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="Organism code")
    parser.add_argument('--ratios', required=True)
    parser.add_argument('--targetdir', required=True)
    parser.add_argument('--numruns', type=int, default=4)
    parser.add_argument('--ncbi_code', default="")
    parser.add_argument('--mincols', type=int, default=8)
    parser.add_argument('--num_cores', type=int, default=1)
    parser.add_argument('--user', default=None)
    parser.add_argument('--csh', action='store_true')
    parser.add_argument('--blocks', default=None)
    parser.add_argument('--inclusion', default=None)
    parser.add_argument('--exclusion', default=None)
    parser.add_argument('--pipeline', default=None)
    parser.add_argument('--setenrich', default=None)
    parser.add_argument('--setenrich_files', default=None)
    args = parser.parse_args()

    if args.csh:
        header = QSUB_TEMPLATE_HEADER_CSH
        template = QSUB_TEMPLATE_CSH
    else:
        header = QSUB_TEMPLATE_HEADER
        template = QSUB_TEMPLATE

    if not os.path.exists(args.targetdir):
      os.makedirs(args.targetdir)

    # write ratios files
    print "Choosing ensemble conditions"
    if args.blocks is None:
      # if inclusion/exlcusion blocks are not defined, simply choose at random
      #
      # I don't think this is right! - needs to be fixed
      #
      dm.prepare_ensemble_matrix( args.ratios, args.targetdir, args.numfiles,
                                 args.mincols )
    else:
      cols = ensemblePicker( ratios = args.ratios, blocks = args.blocks, inclusion = args.inclusion, exclusion = args.exclusion, nruns = args.numruns, ratios_file = args.targetdir )
      cols.pickCols_all()

    # write config files
    config_params = { }

    print "Writing ensemble config files"
    for i in range( 1, args.numruns+1 ):
      if args.setenrich is not None:
        # combine all possible setenrichment modes
        sets = args.setenrich.split(",")
        set_files = args.setenrich_files.split(",")
        
        set_file_ref = {}
        for j in range( 0, len( sets ) ):
          set_file_ref[ sets[ j ] ] = set_files[ j ]

        # don't know how to do it without se
        setcombinations = [ None ]
 
        for j in range( 1, len( sets ) + 1 ):
          setcombinations = setcombinations + [ ( "," ).join( x ) for x in itertools.combinations( sets, j ) ]

        # choose a set enrichment mode
        set_choice = random.sample(setcombinations,1)[0]
        if set_choice is None:
          ini = cMonkeyIniGen( dict( config_params.items() + [ ( "random_seed", i ) ] ) )
          ini.writeIni( os.path.join( args.targetdir, "config-%03d.ini" % i ) )
        else:
          set_choice_files = (",").join( [ set_file_ref[ x ] for x in set_choice.split( "," ) ] )
          ini = cMonkeyIniGen( dict( config_params.items() + [ ( "random_seed", i ), ( "set_types", set_choice ), ( "set_file", set_choice_files ), ("pipeline_file", args.pipeline ) ] ) )
          ini.writeIni( os.path.join( args.targetdir, "config-%03d.ini" % i ) )
      else:
        ini = cMonkeyIniGen( dict( config_params.items() + [ ( "random_seed", i ) ] ) )
        ini.writeIni( os.path.join( args.targetdir, "config-%03d.ini" % i ) )

    with open(os.path.join(args.targetdir, "%s.sh" % args.organism), 'w') as outfile:
        if args.user is not None:
            login = args.user
        else:
            login = os.getlogin()

        outfile.write(header)
        outfile.write(template % (login, args.numruns, login,
                                  args.num_cores,
                                  args.organism,
                                  os.path.join(args.targetdir, "ratios-$BATCHNUM.tsv"),
                                  os.path.join(args.targetdir, "config-$BATCHNUM.ini"),
                                  "%s-out-$BATCHNUM" % (args.organism),
                                  "%s-out-$BATCHNUM" % (args.organism) ) )
