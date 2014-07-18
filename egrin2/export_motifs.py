import sqlite3
import os

MEME_FILE_HEADER = """MEME version 3.0

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
A %.3f C %.3f G %.3f T %.3f
"""

def export_run_motifs_to_meme(dbfile, targetdir, basename):
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    cursor2 = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    iteration = cursor.fetchone()[0]
    filename = '%s.meme' % basename

    a_perc = 0.284
    c_perc = 0.216
    g_perc = 0.216
    t_perc = 0.284

    with open(os.path.join(targetdir, filename), 'w') as outfile:
        outfile.write(MEME_FILE_HEADER % (a_perc, c_perc, g_perc, t_perc))
        cursor.execute('select mi.rowid, cluster, motif_num, evalue, count(mms.pvalue) as num_sites from motif_infos mi join meme_motif_sites mms on mi.rowid = mms.motif_info_id where mi.iteration=? group by mi.rowid', [iteration])
        for rowid, cluster, motif_num, evalue, num_sites in cursor.fetchall():
            motif_name = '%s_%03d_%02d' % (basename, cluster, motif_num)
            outfile.write('\nMOTIF %s\n' % motif_name)
            outfile.write('BL   %s width=0 seqs=0\n' % motif_name)

            cursor2.execute('select a,c,g,t from motif_pssm_rows where motif_info_id=? order by row', [rowid])
            pssm_rows = [(a, c, g, t) for a, c, g, t in cursor2.fetchall()]
            outfile.write('letter-probability matrix: alenght= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
            for a, c, g, t in pssm_rows:
                outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (a, c, g, t))

    cursor.close()
    cursor2.close()
    conn.close()
