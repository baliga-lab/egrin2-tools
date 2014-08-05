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
            motif_name = '%s_%03d_%02d' % (basename, cluster, motif_num)
            outfile.write('\nMOTIF %s\n' % motif_name)
            outfile.write('BL   MOTIF %s width=0 seqs=0\n' % motif_name)

            cursor2.execute('select a,c,g,t from motif_pssm_rows where motif_info_id=? order by row', [rowid])
            pssm_rows = [(a, c, g, t) for a, c, g, t in cursor2.fetchall()]
            outfile.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
            for a, c, g, t in pssm_rows:
                outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (a, c, g, t))

    cursor.close()
    cursor2.close()
    conn.close()
