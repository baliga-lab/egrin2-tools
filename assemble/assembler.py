#!/usr/bin/env python

"""Do all steps to assemble a cMonkey2 ensemble.

Example:

python assembler.py --organism eco --ratios /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz --targetdir /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ --ncbi_code 511145 --ensembledir /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/ --col_annot /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0

python assembler.py --organism eco --ratios ./ratios_eco_m3d.tsv.gz --targetdir ./ --ncbi_code 511145 --ensembledir ./eco-ens-m3d/ --col_annot ./E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0

python assembler.py --organism mtu --ratios ./20141130.MTB.all.ratios.csv.gz --targetdir ./ --ncbi_code 83332 --ensembledir ./  --n_resamples 1000
"""
import argparse
import os, sys, stat
import itertools
import logging
import datetime
import pymongo

import assemble.sql2mongoDB as rdb
from assemble.makeCorems import CoremMaker


QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh

"""

QSUB_TEMPLATE_CSH = """#$ -S /bin/csh
#$ -N %(name)s
#$ -o 'out_messages.txt'
#$ -e 'out_messages.txt'
#$ -m be
#$ -q baliga
#$ -P Bal_%(user)s
#$ -l hostname="baliga2|baliga3"
#$ -M %(user)s@systemsbiology.org
#$ -cwd
#$ -pe serial 1
#$ -l mem_free=4G

python resample.py --host %(host)s --db %(db)s --n_rows %(n_rows)i --n_resamples %(n_resamples)i --port %(port)i

"""

RUN_ALL_TEMPLATE = """#!/bin/bash
for FILE in %s*.sh; do
        qsub $FILE
done
"""

RUN_INFO_TEMPLATE = """[ General ensemble info ]

assembly_date_utc: %(date)s
compiled_by: %(user)s
organism: %(organism)s
ncbi_code: %(ncbi_code)s
host: %(host)s
port: %(port)s
db: %(db)s

[ Corem parameters ]

backbone_pval: %(backbone_pval)s
link_comm_score: %(link_comm_score)s
link_comm_increment: %(link_comm_increment)s
link_comm_density_score: %(link_comm_density_score)s
corem_size_threshold: %(corem_size_threshold)s

[ Condition resamples ]

n_resamples: %(n_resamples)s
"""

DESCRIPTION = """assemble.py - prepare cluster runs"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None # "assembler.log"


def __generate_resample_sge_scripts(db, organism, targetdir, n_resamples, user, host, port):
    logging.debug("connecting to db: '%s'", db.name)
    corem_sizes = list(set([len(i["rows"] ) for i in db["corem"].find({}, {"rows": 1})]))
    logging.debug("corem sizes: '%s'", str(corem_sizes))

    if not os.path.isdir(os.path.abspath(os.path.join(targetdir, "qsub"))):
        os.makedirs(os.path.abspath(os.path.join(targetdir, "qsub")))

    for corem_size in corem_sizes:
        name = organism + "_r_" +  str(corem_size)
        with open(os.path.join(os.path.abspath(os.path.join(targetdir, "qsub")), "%s.sh" % name), 'w') as outfile:
            argss = {
                "user": user,
                "n_resamples": n_resamples,
                "name": name + ".sh",
                "host": host,
                "db": db.name,
                "n_rows": corem_size,
                "port": port
            }
            outfile.write(QSUB_TEMPLATE_HEADER_CSH)
            outfile.write(QSUB_TEMPLATE_CSH % argss)

    with open(os.path.join(os.path.abspath(os.path.join(targetdir, "qsub")), "resample.sh"), 'w') as outfile:
        outfile.write(RUN_ALL_TEMPLATE % (organism + "_r_") )

    logging.info("""Output Qsub scripts to %s.
Transfer these documents to the cluster. Run 'resample.sh' with resample.py in your working directory to compute all resamples.
Once this is done, return here to finish processing corems.""", os.path.abspath(os.path.join(targetdir, "qsub")))
    logging.info("Done.")

    # TODO: Split into 2 parts, because waiting makes this process interactive
    """
    ready = None

    while ready != "Done":
        ready = raw_input("Please type: 'Done' to continue\n")
        """

def __finish_and_dump(corems, resultdb, targetdir):
    corems.finishCorems()
    outfile =  resultdb.prefix + str(datetime.datetime.utcnow()).split(" ")[0] + ".mongodump"
    logging.info("Writing EGRIN2 MongoDB to %s", os.path.join(resultdb.targetdir, outfile))
    add_files = " "
    info = os.path.abspath(os.path.join(targetdir, "ensemble.info"))

    if os.path.isfile(info):
        add_files = add_files + info

    pdf = os.path.join(corems.out_dir, "density_stats.pdf")
    if os.path.isfile(pdf):
        add_files = add_files + pdf
    resultdb.mongo_dump(resultdb.dbname, outfile, add_files=add_files.strip())

    logging.info("Done.")


if __name__ == '__main__':
    import argparse
    import os
    import itertools

    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, type=str, help="3 letter organism code")
    parser.add_argument('--ratios', required=True, help="Path to ratios file. Should be 'raw' (normalized) ratios, not the standardized ratios used by cMonkey")
    parser.add_argument('--targetdir', default='.', required=True, help="Storage path for MongoDB and corem data")
    parser.add_argument('--ncbi_code', required=True, help="NCBI organism code")
    parser.add_argument('--cores', default=3, type=int, help="Number local cores to use for corem C++ scripts")
    parser.add_argument('--ensembledir', default='.', help="Path to ensemble runs. Default: cwd")
    parser.add_argument('--col_annot', default=None, help="Tab-delimited file with experiment annotations")
    parser.add_argument('--host', default="localhost", help="MongoDB host. Default 'localhost'")
    parser.add_argument('--port', default=27017, help="MongoDB port", type=int)
    parser.add_argument('--prefix', default=None, help="Ensemble run prefix. Default: *organism*-out-")
    parser.add_argument('--row_annot', default=None, help="Optional row (gene) annotation tab-delimited file. If not specified, annotations will be downloaded from MicrobesOnline using --ncbi_code.")
    parser.add_argument('--row_annot_match_col', default=None, help="Name of column in row_annot that matches row names in ratios file.")
    parser.add_argument('--gre2motif', default=None, help="Motif->GRE clustering file")
    parser.add_argument('--db', default=None, help="Optional ensemble MongoDB database name")
    parser.add_argument('--genome_annot', default=None, help="Optional genome annotation file. Automatically downloaded from MicrobesOnline using --ncbi_code")
    parser.add_argument('--backbone_pval', default=0.05, type=float, help="Significance pvalue for gene-gene backbone. Default = 0.05.")
    parser.add_argument('--link_comm_score', default=0, type=int, help="Scoring metric for link communities" )
    parser.add_argument('--link_comm_increment', default=0.1, type=float, help="Height increment for cutting agglomerative clustering of link communities" )
    parser.add_argument('--link_comm_density_score', default=5, type=int,  help="Density score for evaluating link communities")
    parser.add_argument('--corem_size_threshold', default=3, type=int, help="Defines minimum corem size. Default = 3." )
    parser.add_argument('--n_resamples', default=10000, type=int, help="Number resamples to compute for corem condition assignment. Default = 10,000")
    parser.add_argument('--cluster', default=True, help="Run re-samples on cluster? Boolean.")
    parser.add_argument('--finish_only', default=False, help="Finish corems only. In case session gets dropped")
    parser.add_argument('--user', default=os.getlogin(), help="Cluster user name")

    args = parser.parse_args()

    dbname = args.db if args.db is not None else "%s_db" % args.organism
    targetdir = args.targetdir
    out_prefix = '%s-out-' % args.organism if args.prefix is None else args.prefix

    info_d = {
        "date": str(datetime.datetime.utcnow()),
        "user": args.user,
        "organism": args.organism,
        "ncbi_code": args.ncbi_code,
        "host": args.host,
        "port": args.port,
        "db": dbname,
        "backbone_pval": args.backbone_pval,
        "link_comm_score": args.link_comm_score,
        "link_comm_increment": args.link_comm_increment,
        "link_comm_density_score": args.link_comm_density_score,
        "corem_size_threshold": args.corem_size_threshold,
        "n_resamples": args.n_resamples
    }
    client = pymongo.MongoClient(host=args.host, port=args.port)

    if not args.finish_only:
        if dbname in client.database_names():
            logging.warn("WARNING: %s database already exists!!!", dbname)
        else:
            logging.info("Initializing MongoDB database: %s", dbname)

    db = client[dbname]

    if not os.path.exists(targetdir):
        os.mkdir(targetdir)

    with open(os.path.abspath(os.path.join(targetdir, "ensemble.info")), 'w') as outfile:
        outfile.write(RUN_INFO_TEMPLATE % info_d)

    resultdb = rdb.ResultDatabase(args.organism, db, args.ensembledir, out_prefix,
                                  args.ratios, args.gre2motif, args.col_annot, args.ncbi_code,
                                  args.genome_annot, args.row_annot,
                                  args.row_annot_match_col, targetdir, db_run_override=False)

    if args.finish_only:
        corems = CoremMaker(args.organism, db, args.backbone_pval, targetdir,
                            args.cores, args.link_comm_score,
                            args.link_comm_increment,
                            args.link_comm_density_score,
                            args.corem_size_threshold,
                            args.n_resamples)

        __finish_and_dump(corems, resultdb, targetdir)
    else:
        if len(resultdb.db_files) > 0:
            resultdb.compile()  # Merge sql into mongoDB
            corems = CoremMaker(args.organism, db, args.backbone_pval, targetdir,
                                args.cores, args.link_comm_score,
                                args.link_comm_increment,
                                args.link_comm_density_score,
                                args.corem_size_threshold,
                                args.n_resamples)
            corems.make_corems()

            if args.cluster:
                __generate_resample_sge_scripts(db, args.organism, targetdir, args.n_resamples, args.user,
                                                args.host, args.port)
            else:
                logging.error("""Non-cluster setup currently not supported. Consider running resamples on a cluster.
This will dramatically speed up this step.""")

            # __finish_and_dump(corems, resultdb, targetdir)
        else:
            logging.error("Could not locate cMonkey2 result files. Please specify --ensembledir")
