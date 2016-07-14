#!/usr/bin/env python3
"""Do all steps to assemble a cMonkey2 ensemble.

Example:

python3 assembler.py --organism eco --ratios /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ratios_eco_m3d.tsv.gz --targetdir /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/ --ncbi_code 511145 --ensembledir /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/eco-ens-m3d/ --col_annot /Users/abrooks/Desktop/Active/Eco_ensemble_python_m3d/E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0

python3 assembler.py --organism eco --ratios ./ratios_eco_m3d.tsv.gz --targetdir ./ --ncbi_code 511145 --ensembledir ./eco-ens-m3d/ --col_annot ./E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0

python3 assembler.py --organism mtu --ratios ./20141130.MTB.all.ratios.csv.gz --targetdir ./ --ncbi_code 83332 --ensembledir ./  --n_resamples 1000


Using the SQLIte engine

PYTHONPATH=. assemble/assembler.py --organism mtu --ratios mtb_files/20141130.MTB.all.ratios.csv --targetdir assemble-py3-test --dbengine sqlite --targetdb mtu-assemble-test.db ass1_input/mtu-out-001/cmonkey_run.db ass1_input/mtu-out-002/cmonkey_run.db ass1_input/mtu-out-004/cmonkey_run.db
"""
import argparse
import os, glob
import itertools
import logging
import datetime
import pymongo
import sqlite3

#import kbase.WorkspaceClient as wsc

import assemble_sqlite as asl
import sql2mongoDB as rdb
from makeCorems import CoremMaker, MongoDB
import resample
import assemble_finish

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


def __generate_resample_sge_scripts(corem_sizes, dbname, organism, targetdir, n_resamples,
                                    user, host, port):
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
                "db": dbname,
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


def make_resample_scripts(args, dbname, targetdir, corem_sizes):
    if args.cluster_arch == 'sge':
        __generate_resample_sge_scripts(corem_sizes, dbname, args.organism,
                                        targetdir, args.n_resamples,
                                        args.sge_user, args.host, args.port)
    else:
        logging.error("""Non-cluster setup currently not supported. Consider running resamples on a cluster.
This will dramatically speed up this step.""")


def merge_sqlite(args):
    if len(args.result_dbs) == 0:
        prefix = '%s-out-' % args.organism if args.prefix is None else args.prefix
        result_dbs = sorted(glob.glob(os.path.join(args.ensembledir, "%s???/cmonkey_run.db" % prefix)))
    else:
        result_dbs = args.result_dbs
    asl.merge(args, result_dbs)
    return True


def merge_mongodb(args, dbclient):
    out_prefix = '%s-out-' % args.organism if args.prefix is None else args.prefix
    resultdb = rdb.ResultDatabase(args.organism, dbclient.dbclient,
                                  args.ensembledir, out_prefix,
                                  args.ratios, args.gre2motif, args.col_annot,
                                  args.genome_annot, args.row_annot,
                                  args.row_annot_match_col, args.targetdir, db_run_override=False)
    if len(resultdb.db_files) > 0:
        resultdb.compile()  # Merge sql into mongoDB
        return True
    else:
        return False


def merge_runs(args, dbclient, dbname):
    info_d = {
        "date": str(datetime.datetime.utcnow()),
        "user": args.sge_user,
        "organism": args.organism,
        "ncbi_code": 0,  #args.ncbi_code,
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

    with open(os.path.abspath(os.path.join(args.targetdir, "ensemble.info")), 'w') as outfile:
        outfile.write(RUN_INFO_TEMPLATE % info_d)

    if args.dbengine == 'mongodb':
        return merge_mongodb(args, dbclient)
    elif args.dbengine == 'sqlite':
        return merge_sqlite(args)
    else:
        raise Exception('Unsupported database engine: %s' % args.dbengine)


def make_corems(args, dbclient):
    """db: use the db client adapter"""
    corems = CoremMaker(args.organism, dbclient, args.backbone_pval, targetdir,
                        args.cores, args.link_comm_score,
                        args.link_comm_increment,
                        args.link_comm_density_score,
                        args.corem_size_threshold,
                        args.n_resamples)
    corems.make_corems()


def make_dbclient(args, dbname):
    if args.dbengine == 'mongodb':
        # connect to MongoDB and check for the database
        client = pymongo.MongoClient(host=args.host, port=args.port)
        if dbname in client.database_names():
            logging.warn("WARNING: %s database already exists!!!", dbname)
        else:
            logging.info("Initializing MongoDB database: %s", dbname)
        db = client[dbname]
        return MongoDB(db)
    elif args.dbengine == 'sqlite':
        conn = sqlite3.connect(dbname)
        return asl.SqliteDB(conn)
    else:
        raise Exception('unknown dbengine: %s' % args.dbengine)


"""
def store_kb_workspace(conn):
    if ('KB_AUTH_TOKEN' in os.environ and 'TARGET_WS' in os.environ and
        'WS_URL' in os.environ):
        auth_token = os.environ['KB_AUTH_TOKEN']
        target_ws = os.environ['TARGET_WS']
        ws_url = os.environ['WS_URL']
        ws_service = wsc.Workspace(ws_url, token=auth_token)
        ws_datatype = 'TODO-0.1'
        data = None
        result = ws_service.save_object({'type': ws_datatype,
                                     'data': data,
                                         'workspace': target_ws})
"""

if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="3 letter organism code")
    parser.add_argument('--ratios', required=True, help="Path to ratios file. Should be 'raw' (normalized) ratios, not the standardized ratios used by cMonkey")
    parser.add_argument('--targetdir', default='.', help="Storage path for MongoDB and corem data")
    parser.add_argument('--backbone_pval', default=0.05, type=float, help="Significance pvalue for gene-gene backbone. Default = 0.05.")
    parser.add_argument('--cores', default=3, type=int, help="Number local cores to use for corem C++ scripts")
    parser.add_argument('--link_comm_score', default=0, type=int, help="Scoring metric for link communities" )
    parser.add_argument('--link_comm_increment', default=0.1, type=float, help="Height increment for cutting agglomerative clustering of link communities" )
    parser.add_argument('--link_comm_density_score', default=5, type=int,  help="Density score for evaluating link communities")
    parser.add_argument('--corem_size_threshold', default=3, type=int, help="Defines minimum corem size. Default = 3." )
    parser.add_argument('--n_resamples', default=10000, type=int, help="Number resamples to compute for corem condition assignment. Default = 10,000")

    #  options for running on Sun Grid Engine
    parser.add_argument('--cluster_arch', default='sge', help="where to run resampling on")
    parser.add_argument('--sge_user', default=os.environ['LOGNAME'], help="Cluster user name")

    parser.add_argument('--dbengine', default='mongodb', help="mongodb or sqlite")

    # MongoDB specific
    parser.add_argument('--host', default="localhost", help="MongoDB host. Default 'localhost'")
    parser.add_argument('--port', default=27017, help="MongoDB port", type=int)
    parser.add_argument('--targetdb', default=None, help="Optional ensemble MongoDB database name")

    # reading from an ensemble directory using directory pattern
    parser.add_argument('--prefix', default=None, help="Ensemble run prefix. Default: *organism*-out-")
    parser.add_argument('--ensembledir', default='.', help="Path to ensemble runs. Default: cwd")

    # These arguments are optional steps
    parser.add_argument('--col_annot', default=None, help="Tab-delimited file with experiment annotations")
    parser.add_argument('--row_annot', default=None, help="Optional row (gene) annotation tab-delimited file. If not specified, annotations will be downloaded from MicrobesOnline using --ncbi_code.")
    parser.add_argument('--row_annot_match_col', default=None, help="Name of column in row_annot that matches row names in ratios file.")
    parser.add_argument('--gre2motif', default=None, help="Motif->GRE clustering file")
    parser.add_argument('--genome_annot', default=None, help="Optional genome annotation file. Automatically downloaded from MicrobesOnline using --ncbi_code")

    # This is alternatively to ensembledir, used by the sqlite merger:
    # specify the input databases individually
    parser.add_argument('result_dbs', nargs='*')
    args = parser.parse_args()

    dbname = args.targetdb if args.targetdb is not None else "%s_db" % args.organism
    targetdir = args.targetdir
    if not os.path.exists(targetdir):
        os.mkdir(targetdir)
    dbclient = make_dbclient(args, dbname)

    if merge_runs(args, dbclient, dbname):
        make_corems(args, dbclient)
        corem_sizes = dbclient.corem_sizes()
        logging.debug("corem sizes: '%s'", str(corem_sizes))

        # when we run on sqlite, we will automatically run resampling
        # right after
        if args.dbengine == 'sqlite':
            rs_dbclient = resample.SqliteDB(dbclient.conn)
            for corem_size in corem_sizes:
                resample.col_resample_ind(rs_dbclient, n_rows=corem_size,
                                          cols=rs_dbclient.get_col_nums(),
                                          n_resamples=1000,
                                          keepP=0.1)
            dbclient.conn.commit()  # commit everything so we can recover here

            # do the finish step
            assemble_finish.finish_corems(dbclient)
            dbclient.conn.commit()
            #store_kb_workspace(dbclient.conn)
        else:
            make_resample_scripts(args, dbname, targetdir, corem_sizes)
    else:
        logging.error("Could not locate cMonkey2 result files. Please specify --ensembledir")
