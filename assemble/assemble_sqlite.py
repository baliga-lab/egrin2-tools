import os
import sqlite3
import logging
import gzip
import pandas as pd
import requests
from datetime import datetime
import itertools
import json

"""
This assemble module is the sqlite3 based implementation
"""

class SqliteDB:
    """database interface to sqlite"""
    def __init__(self, conn):
        self.conn = conn
        self.conn.execute('create table if not exists row_row (keyrow_id int, subrow_id int, counts int, weight decimal, backbone_pval decimal)')
        self.conn.execute('create table if not exists corems (density decimal, corem_num int, weighted_density decimal)')
        self.conn.execute('create table if not exists corem_rows (corem_id int, row_id int)')
        self.conn.execute('create table if not exists corem_cols (corem_id int, col_id int, pval decimal)')
        self.conn.execute('create table if not exists corem_edges (corem_id int, row1_id int, row2_id int)')

    def close(self):
        self.conn.close()

    def get_row_maps(self):
        row2id = {}
        id2row = {}
        cursor = self.conn.cursor()
        try:
            cursor.execute('select rowid,name from rows')
            for rowid, name in cursor.fetchall():
                row2id[name] = rowid
                id2row[rowid] = name
            return row2id, id2row
        finally:
            cursor.close()

    def get_column_maps(self):
        col2id = {}
        id2col = {}
        cursor = self.conn.cursor()
        try:
            cursor.execute('select rowid,name from columns')
            for rowid, name in cursor.fetchall():
                col2id[name] = rowid
                id2col[rowid] = name
            return col2id, id2col
        finally:
            cursor.close()

    def num_row_co_occurence(self, rowname, row2id, id2row):
        """for a given row, return the number of co-occurences with genes in all
        biclusters. The result is a Series using the gene names as indexes and
        counts as values"""
        row_pk = row2id[rowname]
        cursor = self.conn.cursor()
        try:
            cursor.execute('select row_id, count(cluster_id) from bicluster_rows where cluster_id in (select cluster_id from bicluster_rows where row_id=?) group by row_id', [row_pk])
            result = {id2row[row_pk]: count for row_pk, count in cursor.fetchall()}
            return pd.Series(result)
        finally:
            cursor.close()

    def drop_row_rows(self):
        self.conn.execute('delete from row_row')

    def update_row_row(self, keyrow_pk, subrow_pk, data_counts_norm, backbone_pval):
        self.conn.execute('update row_row set weight=?, backbone_pval=? where keyrow_id=? and subrow_id=?', [data_counts_norm, backbone_pval, keyrow_pk, subrow_pk])

    def insert_row_row(self, rowrow_docs):
        for doc in rowrow_docs:
            keyrow_id, subrow_id = map(int, doc['row_ids'])
            counts = int(doc['counts'])
            weight = float(doc['weight'])
            backbone_pval = float(doc['backbone_pval'])
            self.conn.execute('insert into row_row (keyrow_id,subrow_id,counts,weight,backbone_pval) values (?,?,?,?,?)',
                              [keyrow_id, subrow_id, counts, weight, backbone_pval])

    def insert_corem(self, corem_docs):
        cursor = self.conn.cursor()
        try:
            for doc in corem_docs:
                cursor.execute('insert into corems (density,corem_num,weighted_density) values (?,?,?)',
                               [float(doc['density']), int(doc['corem_id']), float(doc['weighted_density'])])
                corem_id = cursor.lastrowid
                for row_id in doc['rows']:
                    cursor.execute('insert into corem_rows (corem_id,row_id) values (?,?)',
                                   [corem_id, row_id])
                for col_id in doc['cols']:
                    cursor.execute('insert into corem_cols (corem_id,col_id) values (?,?)',
                                   [corem_id, col_id])
                for edge_str in doc['edges']:
                    row1, row2 = edge_str.split('-')
                    cursor.execute('insert into corem_edges (corem_id,row1_id,row2_id) values (?,?,?)',
                                   [corem_id, int(row1), int(row2)])
        finally:
            cursor.close()

    def corem_sizes(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute('select distinct count(row_id) as corem_size from corem_rows group by corem_id order by corem_size')
            return [row[0] for row in cursor.fetchall()]
        finally:
            cursor.close()

    def get_cond_ids(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute('select rowid from columns')
            return [row[0] for row in cursor.fetchall()]
        finally:
            cursor.close()

    def get_corems(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute('select corem_id, row_id from corem_rows')
            corem_rows = [(corem_id, row_id) for corem_id, row_id in cursor.fetchall()]
            result = [{'_id': corem_id, 'corem_id': corem_id,
                       'rows': list(map(lambda x: x[1], row_ids))}
                      for corem_id, row_ids in itertools.groupby(corem_rows, lambda x: x[0])]
            return result
        finally:
            cursor.close()

    def no_col_resamples(self, col_id, nrows, nresamples):
        cursor = self.conn.cursor()
        try:
            cursor.execute("select count(*) from col_resamples where col_id=? and nrows=? and nresamples >= ?",
                           [col_id, nrows, nresamples])
            return cursor.fetchone()[0] == 0
        finally:
            cursor.close()

    def find_gene_expressions(self, row_pks, column_pks):
        cursor = self.conn.cursor()
        try:
            row_in_list = '(%s)' % ','.join(map(str, row_pks))
            col_in_list = '(%s)' % ','.join(map(str, column_pks))
            query = 'select col_id,value,std_value from expr_values where row_id in %s and col_id in %s' % (row_in_list, col_in_list)
            cursor.execute(query)
            return pd.DataFrame([{'col_id':  col_id, 'raw_expression': value, 'standardized_expression': std_value}
                                 for col_id, value, std_value in cursor.fetchall()])
        finally:
            cursor.close()

    def find_col_resamples(self, nrows, col_pks):
        cursor = self.conn.cursor()
        try:
            in_list = '(%s)' % ','.join(map(str, col_pks))
            query = 'select col_id,nresamples,lowest_raw_exps,lowest_std_exps from col_resamples where nrows=? and col_id in ' + in_list
            cursor.execute(query, [nrows])
            result = []
            for col_id, nresamples, lowest_raw, lowest_std in cursor.fetchall():
                result.append({'col_id': col_id, 'n_rows': nrows,
                               'resamples': nresamples,
                               'lowest_raw': json.loads(lowest_raw),
                               'lowest_standardized': json.loads(lowest_std)})
            return pd.DataFrame(result)
        finally:
            cursor.close()

    def update_corem(self, corem, new_cols):
        corem_pk = corem['corem_id']
        self.conn.execute('delete from corem_cols where corem_id=?', [corem_pk])
        for col in new_cols:
            self.conn.execute('insert into corem_cols (corem_id,col_id,pval) values (?,?,?)',
                              [corem_pk, int(col['col_id']), col['pval']])



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
    """reads the specified gene expression and returns both
    raw matrix and standardized matrix"""
    if path.endswith('gz'):
        ratios = pd.read_csv(gzip.open(path, 'rb'), index_col=0, sep="\t")
        if ratios.shape[1] == 0:  # attempt using comma as a delimiter if tab failed
            ratios = pd.read_csv(gzip.open(path, 'rb'), index_col=0, sep=",")
    else:
        ratios = pd.read_csv(path, index_col=0, sep="\t")
        if ratios.shape[1] == 0:  # attempt using comma as a delimiter if tab failed
            ratios = pd.read_csv(path, index_col=0, sep=",")

    if ratios.shape[1] == 0:
        raise Exception("Cannot read ratios file. Check delimiter. Should be '\t' or ',' ")
    return ratios, standardize_ratios(ratios)


def create_tables(conn):
    """Create tables in the result database"""
    conn.execute('create table if not exists rows (name text)')
    conn.execute('create table if not exists row_annotations (name text)')
    conn.execute('create table if not exists row_annotation_values (row_id int, annot_id int, value text)')

    conn.execute('create table if not exists columns (name text)')
    conn.execute('create table if not exists col_annotations (name text)')
    conn.execute('create table if not exists col_annotation_values (col_id int, annot_id int, value text)')

    # holds both original and standardized values
    conn.execute('create table if not exists expr_values (row_id int,col_id,value decimal,std_value decimal)')

    # information about individual ensemble runs
    conn.execute('create table if not exists ensemble_runs (date_added timestamp,start_time timestamp,finish_time timestamp,num_iterations int,organism text,species text,num_rows int,num_columns,num_clusters int,git_sha text)')
    conn.execute('create table if not exists ensemble_run_rows (run_id int, row_id int)')
    conn.execute('create table if not exists ensemble_run_cols (run_id int, col_id int)')

    # bicluster information
    conn.execute('create table if not exists biclusters (run_id int, cluster_num int, residual decimal)')
    conn.execute('create table if not exists bicluster_rows (cluster_id int, row_id int)')
    conn.execute('create table if not exists bicluster_cols (cluster_id int, col_id int)')

    conn.execute('create table if not exists motif_infos (cluster_id int, seqtype text, motif_num int, evalue decimal)')
    conn.execute('create table if not exists motif_pssm_rows (motif_info_id int, row int, a decimal, c decimal, g decimal, t decimal)')
    conn.execute('create table if not exists meme_motif_sites (motif_info_id int, seq_name text, reverse boolean, start int, pvalue decimal)')

    # indexes
    conn.execute('create index if not exists rows_idx on rows (name)')
    conn.execute('create index if not exists row_annotations_idx on row_annotations (name)')

    conn.execute('create index if not exists cols_idx on columns (name)')
    conn.execute('create index if not exists col_annotations_idx on col_annotations (name)')

def create_expr_indexes(conn):
    conn.execute('create index if not exists expr_values_idx on expr_values (row_id,col_id)')


def annotate_microbes_online(conn, row2id, ncbi_code):
    """add all microbes online columns as an attribute to the row"""
    resp = requests.get("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=tab" % ncbi_code)
    content = [line.split('\t') for line in resp.text.split('\n')]
    titles = content[0]
    sysname_col = titles.index('sysName')
    print("sys name column at index: ", sysname_col)

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


def store_ratios(conn, raw_ratios, std_ratios, row2id, col2id):
    """Store gene expressions"""
    logging.info("Storing gene expressions...")
    num_rows = 0

    # speed up access by storing frequent references in variables
    raw_ratios_colnames = raw_ratios.columns.values
    raw_ratio_vals = raw_ratios.values
    std_ratio_vals = std_ratios.values

    for rowidx, rowname in enumerate(raw_ratios.index.values):
        if num_rows % 200 == 0:
            logging.info("%.2f percent done (%d rows)",
                         round((float(num_rows) / raw_ratios.shape[0]) * 100, 1), num_rows)
        row_pk = row2id[rowname]
        to_insert = []
        for colidx, colname in enumerate(raw_ratios_colnames):
            col_pk = col2id[colname]
            raw_value = raw_ratio_vals[rowidx, colidx]
            std_value = std_ratio_vals[rowidx, colidx]
            to_insert.append((row_pk, col_pk, raw_value, std_value))
        conn.executemany('insert into expr_values (row_id,col_id,value,std_value) values (?,?,?,?)',
                         to_insert)
        num_rows += 1

    conn.commit()
    logging.info("done.")


def store_run_info(conn, src_conn, row2id, col2id):
    """Stores the information about an individual ensemble run in the database"""
    logging.info("Store individual run information...")
    src_cur = src_conn.cursor()
    cursor = conn.cursor()
    try:
        src_cur.execute('select start_time,finish_time,num_iterations,organism,species,num_rows,num_columns,num_clusters,git_sha from run_infos')
        run_info = src_cur.fetchone()
        cursor.execute('insert into ensemble_runs (date_added,start_time,finish_time,num_iterations,organism,species,num_rows,num_columns,num_clusters,git_sha) values (?,?,?,?,?,?,?,?,?,?)',
                        [datetime.now(), run_info[0], run_info[1], run_info[2],
                        run_info[3], run_info[4], run_info[5], run_info[6],
                        run_info[7], run_info[8]])
        run_id = cursor.lastrowid
        src_cur.execute('select name from row_names')
        row_names = [row[0] for row in src_cur.fetchall()]
        src_cur.execute('select name from column_names')
        col_names = [row[0] for row in src_cur.fetchall()]
        for rowname in row_names:
            cursor.execute('insert into ensemble_run_rows (run_id,row_id) values (?,?)',
                           [run_id, row2id[rowname]])
        for colname in col_names:
            cursor.execute('insert into ensemble_run_cols (run_id,col_id) values (?,?)',
                           [run_id, col2id[colname]])

        conn.commit()
        return run_id
    finally:
        cursor.close()
        src_cur.close()


def store_biclusters(conn, src_conn, run_id, row2id, col2id):
    """copy bicluster information for the specified run"""
    logging.info("copying biclusters...")
    src_cursor = src_conn.cursor()
    src_cursor2 = src_conn.cursor()
    cursor = conn.cursor()
    cluster2id = {}
    try:
        src_cursor.execute('select max(iteration) from cluster_stats')
        last_iter = src_cursor.fetchone()[0]
        src_cursor.execute('select cluster,residual from cluster_stats where iteration=?',
                           [last_iter])
        for cluster, residual in src_cursor.fetchall():
            src_cursor2.execute('select name from row_members rm join row_names rn on rm.order_num=rn.order_num where rm.iteration=? and rm.cluster=?', [last_iter, cluster])
            rownames = [row[0] for row in src_cursor2.fetchall()]
            src_cursor2.execute('select name from column_members cm join column_names cn on cm.order_num=cn.order_num where cm.iteration=? and cm.cluster=?', [last_iter, cluster])
            colnames = [row[0] for row in src_cursor2.fetchall()]

            cursor.execute('insert into biclusters (run_id,cluster_num,residual) values (?,?,?)',
                           [run_id, cluster, residual])
            cluster_id = cursor.lastrowid
            cluster2id[cluster] = cluster_id

            for rowname in rownames:
                cursor.execute('insert into bicluster_rows (cluster_id,row_id) values (?,?)',
                               [cluster_id, row2id[rowname]])
            for colname in colnames:
                cursor.execute('insert into bicluster_cols (cluster_id,col_id) values (?,?)',
                               [cluster_id, col2id[colname]])

        conn.commit()
        return cluster2id
    finally:
        cursor.close()
        src_cursor.close()
        src_cursor2.close()


def store_motifs(conn, src_conn, cluster2id):
    src_cursor = src_conn.cursor()
    src_cursor2 = src_conn.cursor()
    cursor = conn.cursor()
    try:
        src_cursor.execute('select max(iteration) from cluster_stats')
        last_iter = src_cursor.fetchone()[0]
        src_cursor.execute('select rowid,cluster,seqtype,motif_num,evalue from motif_infos where iteration=?', [last_iter])
        for motif_info_id, cluster, seqtype, motif_num, evalue in src_cursor.fetchall():
            cursor.execute('insert into motif_infos (cluster_id,seqtype,motif_num,evalue) values (?,?,?,?)', [cluster2id[cluster], seqtype, motif_num, evalue])
            motif_id = cursor.lastrowid

            # PSSMs
            src_cursor2.execute('select row,a,c,g,t from motif_pssm_rows where motif_info_id=? order by row', [motif_info_id])
            for row, a, c, g, t in src_cursor2.fetchall():
                cursor.execute('insert into motif_pssm_rows (motif_info_id,row,a,c,g,t) values (?,?,?,?,?,?)', [motif_id, row, a, c, g, t])

            # Sites
            src_cursor2.execute('select seq_name,reverse,start,pvalue from meme_motif_sites where motif_info_id=?', [motif_info_id])
            for seqname, reverse, start, pvalue in src_cursor2.fetchall():
                cursor.execute('insert into meme_motif_sites (motif_info_id,seq_name,reverse,start,pvalue) values (?,?,?,?,?)', [motif_id, seqname, reverse, start, pvalue])

        conn.commit()
    finally:
        src_cursor.close()
        src_cursor2.close()
        cursor.close()


def merge(dbclient, args, result_dbs):
    conn = sqlite3.connect(args.targetdb, 15, isolation_level='DEFERRED')
    create_tables(conn)
    cmonkey_dbs = list(filter(is_valid_db, result_dbs))
    if len(cmonkey_dbs) > 0:
        ncbi_code = extract_ncbi_code(cmonkey_dbs[0])
        print("NCBI code: ", ncbi_code)
        raw_ratios, std_ratios = read_ratios(args.ratios)
        row2id = db_insert_rows(conn, raw_ratios.index.values)
        col2id = db_insert_cols(conn, raw_ratios.columns.values)
        annotate_microbes_online(conn, row2id, ncbi_code)
        conn.commit()  # safe point
        store_ratios(conn, raw_ratios, std_ratios, row2id, col2id)
        conn.commit()  # safe point
        create_expr_indexes(conn)

        for cmonkey_db in cmonkey_dbs:
            src_conn = sqlite3.connect(cmonkey_db)
            try:
                run_id = store_run_info(conn, src_conn, row2id, col2id)
                cluster2id = store_biclusters(conn, src_conn, run_id, row2id, col2id)
                store_motifs(conn, src_conn, cluster2id)
            finally:
                src_conn.close()
    else:
        raise Exception('no input databases provided !!!')
