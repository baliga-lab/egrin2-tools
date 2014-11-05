import numpy as np
import os
import pandas as pd
import weblogolib as wl
import cStringIO
import sqlite3

from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def system( cmd ):
    print cmd
    tmp = os.popen( cmd ).read()
    return tmp

def stop():
    raise 'STOPPING ON PURPOSE!'

def get_motif_pssm(dbfile, cluster_num, motif_num):
    """export the specified motif to a pandas dataframe
    Parameters:
    - dbfile: path to the result database file
    - cluster_num: bicluster number
    - motif_num: motif number
    """
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    iteration = cursor.fetchone()[0]

    query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
    params = [iteration, cluster_num, motif_num]
    cursor.execute(query, params)
    rowid = cursor.fetchone()[0]

    query = 'select a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?'
    params = [iteration, rowid]
    pssm = pd.read_sql( query, conn, params=params )
    return pssm

def get_biop_motif(dbfile, cluster_num, motif_num, option='sites'):
    ##import egrin2.export_motifs as em
    """export the specified motif to a biopython motif object
    Parameters:
    - dbfile: path to the result database file
    - cluster_num: bicluster number
    - motif_num: motif number
    - option of how to translate - sites: jaspar 'sites' file; pfm: jaspar 'pfm' file
    """
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    iteration = cursor.fetchone()[0]

    query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
    params = [iteration, cluster_num, motif_num]
    cursor.execute(query, params)
    rowid = cursor.fetchone()[0]
    mot_info = pd.read_sql('select * from motif_infos where rowid=?', conn, params=[rowid])
    mot_sites = pd.read_sql('select * from meme_motif_sites where motif_info_id=?', conn, params=[rowid])
    
    output = cStringIO.StringIO()
    ## ONE WAY TO TRY -- but Bio.motifs cant parse the incomplete MEME file
    ##output.write(em.MEME_FILE_HEADER % (0.25, 0.25, 0.25, 0.25))
    ##em.write_pssm(output, cursor, os.path.dirname(dbfile), cluster_num, rowid,
    ##              motif_num, mot_info['evalue'][0], 10)
    ##output.seek(0)
    ##mot = motifs.read( output, 'meme' )
    
    ## Second way - create a jaspar 'pfm' file from the pssm
    if option == 'pfm':
        query = 'select a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?'
        params = [iteration, rowid]
        pssm = pd.read_sql( query, conn, params=params )
        counts = np.round( pssm * mot_sites.shape[0] ).transpose()
        counts.to_string(output, header=False, index=False )
        output.seek(0)
        mot = motifs.read( output, 'pfm' )

    ## Third way - create a jaspar 'sites' file
    elif option == 'sites':
        seqs = {}
        for i in xrange(mot_sites.shape[0]):
            name = mot_sites.ix[i].seq_name
            flank_left = mot_sites.ix[i].flank_left
            flank_left = Seq(flank_left if flank_left is not None else "", IUPAC.IUPACAmbiguousDNA()).lower()
            seq = Seq(mot_sites.ix[i].seq, IUPAC.IUPACAmbiguousDNA())
            flank_right = mot_sites.ix[i].flank_right
            flank_right = Seq(flank_right if flank_right is not None else "", IUPAC.IUPACAmbiguousDNA()).lower()
            full_seq = flank_left + seq + flank_right
            bs = SeqRecord( full_seq, id=name )
            seqs[i] = bs
        
        SeqIO.write(seqs.values(), output, 'fasta')
        output.seek(0)
        mot = motifs.read( output, 'sites' )

    output.close()
    ## mot.weblogo('test.png')
    return mot

## Note Bio.motifs.weblogo() uses the weblogo server (slow? requires connection.)
def plot_pssm( mot ):
    kwargs = dict(color_scheme='classic')
    mot.weblogo('file.png', color_scheme='color_classic') ## note, can use format='PDF'
    img = mpimg.imread('file.png')
    imgplot = plt.imshow( img )

## This uses weblogolib package to create files directly (installed as weblogo via pip)
## https://code.google.com/p/weblogo/
def plot_motif( dbfile, cluster_num, motif_num ):
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    iteration = cursor.fetchone()[0]

    query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
    params = [iteration, cluster_num, motif_num]
    cursor.execute(query, params)
    rowid = cursor.fetchone()[0]
    mot_info = pd.read_sql('select * from motif_infos where rowid=?', conn, params=[rowid])
    mot_sites = pd.read_sql('select * from meme_motif_sites where motif_info_id=?', conn, params=[rowid])

    ldata = wl.LogoData.from_seqs(wl.SeqList(mot_sites.seq.values.tolist(), wl.unambiguous_dna_alphabet))
    options = wl.LogoOptions()
    options.fineprint = os.path.dirname(dbfile) + ' %03d %03d' % ( cluster_num, motif_num )
    format = wl.LogoFormat(ldata, options) 
    format.color_scheme = wl.classic
    format.resolution = 150
    tmp = wl.png_formatter( ldata, format )

    output = cStringIO.StringIO(tmp)
    img = mpimg.imread(output)
    imgplot = plt.imshow( img )
