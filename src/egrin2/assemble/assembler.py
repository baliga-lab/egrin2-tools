#!/usr/bin/env python3
"""Do all steps to assemble a cMonkey2 ensemble.
"""
import argparse
import os, glob
import itertools
import logging
import datetime
import sqlite3

import egrin2.assemble.assemble_sqlite as asl
from egrin2.assemble.makeCorems import CoremMaker
import egrin2.assemble.resample as resample
import egrin2.assemble.assemble_finish as assemble_finish
import egrin2.assemble.motif_clusters as motif_clusters


RUN_INFO_TEMPLATE = """[ General ensemble info ]

assembly_date_utc: %(date)s
organism: %(organism)s
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
LOG_FILE = None


def merge_sqlite(args, dbclient):
    if len(args.result_dbs) == 0:
        prefix = '%s-out-' % args.organism if args.prefix is None else args.prefix
        result_dbs = sorted(glob.glob(os.path.join(args.ensembledir, "%s???/cmonkey_run.db" % prefix)))
    else:
        print('args.result_dbs: ', args.result_dbs)
        result_dbs = args.result_dbs
    asl.merge(dbclient, args, result_dbs)
    return True


def merge_runs(args, dbclient, dbname):
    info_d = {
        "date": str(datetime.datetime.utcnow()),
        "organism": args.organism,
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

    return merge_sqlite(args, dbclient)


def make_corems(args, dbclient):
    """db: use the db client adapter"""
    corems = CoremMaker(args.organism, dbclient, args.backbone_pval, args.targetdir,
                        args.cores, args.link_comm_score,
                        args.link_comm_increment,
                        args.link_comm_density_score,
                        args.corem_size_threshold,
                        args.n_resamples)
    corems.make_corems()


def make_dbclient(targetdb):
    if targetdb is None:
        raise Exception("no sqlite database file specified !")
    conn = sqlite3.connect(targetdb, 15, isolation_level='DEFERRED')
    return asl.SqliteDB(conn)


def main():
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

    parser.add_argument('--targetdb', required=True, help="Sqlite3 result database name to be created")

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

    targetdir = args.targetdir
    if not os.path.exists(targetdir):
        os.mkdir(targetdir)
    dbclient = make_dbclient(args.targetdb)

    if merge_runs(args, dbclient, args.targetdb):
        make_corems(args, dbclient)
        corem_sizes = dbclient.corem_sizes()
        logging.debug("corem sizes: '%s'", str(corem_sizes))

        # run the whole assemble step
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

        # store GRE information into motif_infos table
        motif2gre = motif_clusters.load_gre_map(args.ensembledir)
        motif_clusters.store_gres(dbclient.conn, motif2gre)  # does commit

        # TODO: insert FIMO (fimo tables fimo and fimo_small)
    else:
        logging.error("Could not locate cMonkey2 result files. Please specify --ensembledir")


if __name__ == '__main__':
    main()
