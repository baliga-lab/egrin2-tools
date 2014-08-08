import sqlite3
import os

"""export_motifs.py - module to support motif extraction from
a cmonkey-python run database
"""

MEME_FILE_HEADER = """MEME version 3.0

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
A %.3f C %.3f G %.3f T %.3f
"""

def write_pssm(outfile, cursor, run_name, cluster, motif_info_id,
               motif_num, evalue, num_sites):
    """writes a single PSSM to the given file"""
    motif_name = '%s_%03d_%02d' % (run_name, cluster, motif_num)
    outfile.write('\nMOTIF %s\n' % motif_name)
    outfile.write('BL   MOTIF %s width=0 seqs=0\n' % motif_name)

    cursor.execute('select a,c,g,t from motif_pssm_rows where motif_info_id=? order by row',
                   [motif_info_id])
    pssm_rows = [(a, c, g, t) for a, c, g, t in cursor.fetchall()]
    outfile.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
    for a, c, g, t in pssm_rows:
        outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (a, c, g, t))


def export_run_motifs_to_meme(dbfile, targetdir, basename, max_residual=None,
                              max_evalue=None):
    """export the specified run's motifs to a file in MEME format.
    Parameters:
    - dbfile: path to the result database file
    - targetdir: target directory for the output file
    - basename: a base name for the output which is used as a prefix
    """
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    cursor2 = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    iteration = cursor.fetchone()[0]
    filename = '%s.meme' % basename

    # these are currently just for legacy runs, cmonkey-python has a table for
    # global background now
    a_perc = 0.284
    c_perc = 0.216
    g_perc = 0.216
    t_perc = 0.284

    num_pssms = 0
    with open(os.path.join(targetdir, filename), 'w') as outfile:
        outfile.write(MEME_FILE_HEADER % (a_perc, c_perc, g_perc, t_perc))
        query = 'select mi.rowid, mi.cluster, motif_num, evalue, count(mms.pvalue) as num_sites from motif_infos mi join meme_motif_sites mms on mi.rowid = mms.motif_info_id join cluster_residuals cr on mi.cluster=cr.cluster and mi.iteration=cr.iteration where mi.iteration=?'
        params = [iteration]
        if max_residual is not None:
            query += ' and cr.residual <= ?'
            params.append(max_residual)
        if max_evalue is not None:
            query += ' and mi.evalue <= ?'
            params.append(max_evalue)
        query += ' group by mi.rowid'
        cursor.execute(query, params)
        for rowid, cluster, motif_num, evalue, num_sites in cursor.fetchall():
            write_pssm(outfile, cursor2, basename, cluster, rowid, motif_num, evalue,
                       num_sites)
            num_pssms += 1

    cursor.close()
    cursor2.close()
    conn.close()
    return num_pssms


def make_meme_file(dbpaths, maxiter, targetdir, gene,
                   max_residual=None, max_evalue=None,
                   a_perc=0.284, c_perc=0.216, g_perc=0.216, t_perc=0.284):
    """Creates a meme file for a specific gene. Returns the number of
    PSSMs that were written"""
    num_written = 0
    targetpath = os.path.join(targetdir, '%s.meme' % gene)
    with open(targetpath, 'w') as outfile:
        outfile.write(MEME_FILE_HEADER % (a_perc, c_perc, g_perc, t_perc))
        # gene -> all runs
        for dbpath in dbpaths:
            conn = sqlite3.connect(dbpath)
            cursor = conn.cursor()
            cursor2 = conn.cursor()
            cursor.execute('select order_num from row_names where name=?', [gene])
            order_num = cursor.fetchone()[0]
            query = """select mi.rowid, mi.cluster, motif_num, evalue,
                       count(mms.pvalue) as num_sites
                       from (select rowid,cluster,motif_num,evalue from motif_infos
                             where iteration=?) mi
                       join (select cluster from row_members where iteration=?
                              and order_num=?) as rm on mi.cluster=rm.cluster
                       join (select cluster, residual from cluster_residuals
                             where iteration=?) cr on rm.cluster=cr.cluster
                       join meme_motif_sites mms on mi.rowid=mms.motif_info_id"""
            params = [maxiter, maxiter, order_num, maxiter]

            # optional residual and evalue filters
            if max_residual is not None or max_evalue is not None:
                query += " where"
            if max_residual is not None:
                query += " residual <= ?"
                params.append(max_residual)
            if max_evalue is not None:
                if max_residual is not None:
                    query += " and "
                query += " evalue <= ?"
                params.append(max_evalue)

            query += " group by mi.rowid"
            cursor.execute(query, params)
            for rowid, cluster, motif_num, evalue, num_sites in cursor.fetchall():
                write_pssm(outfile, cursor2,
                           os.path.basename(os.path.dirname(dbpath)),
                           cluster, rowid, motif_num, evalue, num_sites)
                num_written += 1
            cursor2.close()
            cursor.close()
            conn.close()
    return num_written


def make_meme_files(inpath, prefix, targetdir):
    """create MEME files based on each gene in the ensemble run and writing
    all clusters in all runs that contain the gene"""
    resultdirs = map(lambda s: os.path.join(inpath, s),
                     sorted([entry for entry in os.listdir(inpath)
                             if entry.startswith(prefix) and os.path.isdir(entry)]))
    dbpaths = [os.path.join(resultdir, 'cmonkey_run.db') for resultdir in resultdirs]
    # extract max iteration
    conn = sqlite3.connect(dbpaths[0])
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from row_members')
    max_iteration = cursor.fetchone()[0]
    cursor.execute('select name from row_names order by name')
    genes = [row[0] for row in cursor.fetchall()]
    conn.close()
    for gene in genes:
        print "processing gene '%s'..." % gene,
        num_written = make_meme_file(dbpaths, max_iteration, targetdir, gene)
        print "%d motifs written." % num_written
