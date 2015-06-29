#!/usr/bin/env python

import argparse
import os
import sqlite3
import logging
import gzip
import pandas as pd
import requests

"""
This merger tool outfactors the merging process, and merges into
sqlite instead to reduce dependency
"""

DESCRIPTION = """merger.py - merge cmonkey2 ensemble runs"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = None


def is_valid_db(dbpath):
    conn = sqlite3.connect(dbpath)
    cur = conn.cursor()
    try:
        cur.execute('select finish_time from run_infos')
        row = cur.fetchone()
        return row is not None and row[0] is not None
    except:
        logging.error("incomplete run")
        return False
    finally:
        cur.close()
        conn.close()


def extract_ncbi_code(dbpath):
    conn = sqlite3.connect(dbpath)
    cur = conn.cursor()
    try:
        cur.execute('select ncbi_code from run_infos')
        return cur.fetchone()[0]
    finally:
        cur.close()
        conn.close()


def standardize_ratios(ratios):
    """compute standardized ratios (global). row standardized"""
    ratios_standardized = ratios.copy()
    zscore = lambda x: (x - x.mean()) / x.std()

    for row in ratios.iterrows():
        ratios_standardized.loc[row[0]] = zscore(row[1])

    return ratios_standardized


def read_ratios(path):    
    ratios = pd.read_csv(gzip.open(path, 'rb'), index_col=0, sep="\t")

    if ratios.shape[1] == 0:  # attempt using comma as a delimiter if tab failed
        ratios = pd.read_csv(gzip.open(path, 'rb'), index_col=0, sep=",")

    if ratios.shape[1] == 0:
        raise Exception("Cannot read ratios file. Check delimiter. Should be '\t' or ',' ")
    return standardize_ratios(ratios)


def create_tables(conn):
    """Create tables in the result database"""
    conn.execute('create table rows (name text)')
    conn.execute('create table row_annotations (name text)')
    conn.execute('create table row_annotation_values (row_id int, annot_id int, value text)')

    conn.execute('create table columns (name text)')
    conn.execute('create table col_annotations (name text)')
    conn.execute('create table col_annotation_values (col_id int, annot_id int, value text)')


def annotate_microbes_online(conn, row2id, ncbi_code):
    """add all microbes online columns as an attribute to the row"""
    resp = requests.get("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=tab" % ncbi_code)
    content = [line.split('\t') for line in resp.text.split('\n')]
    titles = content[0]
    sysname_col = titles.index('sysName')
    print "sys name column at index: ", sysname_col
    
    # insert the attributes
    attr2id = {}
    cursor = conn.cursor()
    try:
        for attribute in [title for title in titles if title != 'sysName']:
            cursor.execute('insert into row_annotations (name) values (?)', [attribute])
            attr2id[attribute] = cursor.lastrowid

        for line in content[1:]:
            if len(line) < sysname_col + 1:  # line too short ?
                continue
            sysname = line[sysname_col]
            print "processing ", sysname
            if sysname in row2id:
                row_pk = row2id[sysname]
                for col_idx in range(len(line)):
                    value = line[col_idx].strip()
                    if col_idx != sysname_col and len(value) > 0:
                        attr_id = attr2id[titles[col_idx]]
                        cursor.execute("insert into row_annotation_values (row_id,annot_id,value) values (?,?,?)", [row2id[sysname], attr_id, line[col_idx]])
        conn.commit()
    except:
        conn.rollback()
        raise
    finally:
        cursor.close()    


def db_insert_rows(conn, rows):
    cursor = conn.cursor()
    result = {}
    for row in rows:
        cursor.execute('insert into rows (name) values (?)', [row])
        result[row] = cursor.lastrowid
    conn.commit()
    return result


def db_insert_cols(conn, cols):
    cursor = conn.cursor()
    result = {}
    for col in cols:
        cursor.execute('insert into columns (name) values (?)', [col])
        result[col] = cursor.lastrowid
    conn.commit()
    return result


if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, type=str,
                        help="3 letter organism code")
    parser.add_argument('--ratios', required=True)
    parser.add_argument('--targetdb', required=True)
    parser.add_argument('result_dbs', nargs='*')
    args = parser.parse_args()
    conn = sqlite3.connect(args.targetdb, 15, isolation_level='DEFERRED')
    #conn = sqlite3.connect(args.targetdb)
    try:
        create_tables(conn)
        cmonkey_dbs = filter(is_valid_db, args.result_dbs)
        print cmonkey_dbs
        if len(cmonkey_dbs) > 0:
            ncbi_code = extract_ncbi_code(cmonkey_dbs[0])
            print "NCBI code: ", ncbi_code
            ratios = read_ratios(args.ratios)
            row2id = db_insert_rows(conn, ratios.index.values)
            col2id = db_insert_cols(conn, ratios.columns.values)
            print col2id
            #row2id = { name: i for i, name in enumerate(ratios.index.values) }
            #col2id = { name: i for i, name in enumerate(ratios.columns.values) }
            annotate_microbes_online(conn, row2id, ncbi_code)
    finally:
        conn.close()

    
"""Mongo RowInfo
{
    "_id" : ObjectId("558c770b77ffbc5e8711a251"),
"sysName" : "Rv0020c",
"TIGRRoles" : NaN,
"scaffoldId" : 7022,
"name" : "TB39.8",
"locusId" : 31791,
"egrin2_row_name" : "Rv0020c",
"COGDesc" : NaN,
"stop" : 23861,
"accession" : "NP_214534.1",
"EC" : NaN,
"start" : 25444,
"COGFun" : NaN,
"strand" : "-",
"COG" : NaN,
"TIGRFam" : NaN,
"row_id" : 19,
"ECDesc" : NaN,
"GO" : NaN,
"GI" : 15607162,
"desc" : "hypothetical protein (NCBI)"
}
"""
"""Mongo ColInfo
{
    "_id" : ObjectId("558c770c77ffbc5e8711b198"),
"col_id" : 19,
"egrin2_col_name" : "H37Rv..Youman.media.replicate.4",
"additional_info" : [ ]
}
"""
"""Mongo BiclusterInfo
{
    "_id" : ObjectId("558c781077ffbc5e877ab582"),
"rows" : [  972,  1810,  2124,  1047,  789,  2400,  265,  3092,  3727 ],
"run_id" : 0,
"residual" : 0.5492853502653435,
"cluster" : 19,
"columns" : [ 	982, 	692, 	693, 	494, 	739, 	737, 	738, 	701, 	984, 	985, 	981, 	292, 	967, 	966, 	746, 	947, 	753, 	752, 	291, 	920, 	290, 	919, 	937, 	485, 	487, 	482, 	483, 	484, 	491, 	490, 	289, 	288, 	159, 	1638, 	1640, 	1641, 	1170, 	344, 	496, 	346, 	347, 	349, 	707, 	705, 	706, 	704, 	731, 	1166, 	730, 	729, 	740, 	741, 	479, 	723, 	721, 	722, 	708, 	709, 	710, 	711, 	1168, 	974, 	972, 	973, 	970, 	971, 	969, 	1167, 	960, 	509, 	747, 	744, 	1164, 	745, 	951, 	952, 	953, 	954, 	1639, 	284, 	748, 	749, 	1169, 	156, 	507, 	503, 	506, 	388, 	714, 	713, 	712, 	286, 	287 ]
}
"""
"""Mongo EnsembleInfo
{ "_id" : ObjectId("558c77f877ffbc5e877ab56f"), "run_id" : 4, "added_to_ensemble" : ISODate("2015-06-25T21:51:52.020Z"), "start_time" : "2015-01-23 13:22:07.160374", "run_name" : "mtu-out-008", "cols" : [ 1642, 	413, 	1685, 	1244, ... ], "num_columns" : 184, "species" : "Mycobacterium_tuberculosis_H37Rv", "git_sha" : "$Id: 69b1afe09739513dfe1102407ea04fe0e5c65e9e $", "rows" : [ 0, 1, 2, 3, 4, 5, ...], "finish_time" : "2015-01-24 11:46:47.914951", "num_iterations" : 2000, "num_rows" : 3923, "num_clusters" : 420, "organism" : "mtu" }
"""
"""Mongo Corem
{ "_id" : ObjectId("558c79f277ffbc5e87859004"), "rows" : [  3219,  3461,  3462 ], "density" : 1, "corem_id" : 20, "cols" : [ ], "edges" : [  "3219-3462",  "3219-3461",  "3462-3461" ], "weighted_density" : 0.03213 }
"""
"""Mongo expression
{ "_id" : ObjectId("558c775277ffbc5e8711b871"), "col_id" : 19, "raw_expression" : 0.667, "row_id" : 0, "standardized_expression" : 0.9155509215108877 }
"""
"""Mongo genome: is DNA sequence"""
"""Mongo motif info
{ "_id" : ObjectId("558c783777ffbc5e877ab78e"), "meme_motif_site" : [ 	{ 	"start" : 86, 	"row_id" : 1602, "scaffoldId" : 7022, 	"reverse" : 0, 	"pvalue" : 8.05e-7 }, 	{ 	"start" : 50, 	"row_id" : 511, 	"scaffoldId" : 7022, 	"reverse" : 0, 	"pvalue" : 0.00000189 }, 	{ 	"start" : 63, 	"row_id" : 289, "scaffoldId" : 7022, 	"reverse" : 1, 	"pvalue" : 0.0000271 }, 	{ 	"start" : 78, 	"row_id" : 49, 	"scaffoldId" : 7022, 	"reverse" : 0, 	"pvalue" : 0.0000316 }, 	{ 	"start" : 59, 	"row_id" : 1121, "scaffoldId" : 7022, 	"reverse" : 0, 	"pvalue" : 0.0000356 }, 	{ 	"start" : 87, 	"row_id" : 2172, "scaffoldId" : 7022, 	"reverse" : 1, 	"pvalue" : 0.0000744 }, 	{ 	"start" : 67, 	"row_id" : 2214, "scaffoldId" : 7022, 	"reverse" : 0, 	"pvalue" : 0.000197 } ], "evalue" : 1900, "gre_id" : "NaN", "cluster_id" : ObjectId("558c781077ffbc5e877ab57a"), "motif_num" : 2, "seqtype" : "upstream", "pwm" : [ 	{ 	"a" : 0.857143, 	"c" : 0, 	"t" : 0.142857, 	"g" : 0, 	"row" : 0 }, 	{ 	"a" : 0.714286, 	"c" : 0, 	"t" : 0.285714, 	"g" : 0, 	"row" : 1 }, 	{ 	"a" : 0.571429, 	"c" : 0, "t" : 0.428571, 	"g" : 0, 	"row" : 2 }, 	{ 	"a" : 0, 	"c" : 0.285714, 	"t" : 0.142857, 	"g" : 0.571429, 	"row" : 3 }, 	{ 	"a" : 0, 	"c" : 0, 	"t" : 0, 	"g" : 1, "row" : 4 }, 	{ 	"a" : 0.142857, 	"c" : 0, 	"t" : 0.428571, 	"g" : 0.428571, 	"row" : 5 }, 	{ 	"a" : 0.142857, 	"c" : 0.714286, 	"t" : 0, 	"g" : 0.142857, 	"row" : 6 }, 	{ 	"a" : 0, 	"c" : 0, 	"t" : 0.714286, 	"g" : 0.285714, 	"row" : 7 }, 	{ 	"a" : 1, 	"c" : 0, 	"t" : 0, 	"g" : 0, 	"row" : 8 }, 	{ 	"a" : 0, "c" : 0.857143, 	"t" : 0, 	"g" : 0.142857, 	"row" : 9 } ] }
"""
"""Mongo rowrow
{ "_id" : ObjectId("558c798177ffbc5e877acd5e"), "row_ids" : [  238,  3779 ], "counts" : 2, "backbone_pval" : 0.17728185409972097, "weight" : 0.009950248756218905 }
"""
