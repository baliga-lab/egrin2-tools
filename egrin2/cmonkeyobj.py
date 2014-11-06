import numpy as np
import os
import pandas as pd
import cStringIO
import sqlite3
import bz2
import cPickle as pickle

from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import weblogolib as wl

def read_all_tables( dbfile ):
    """read out all tables in a sqlite3 db file into a dict of pandas dataframes"""
    conn = sqlite3.connect( dbfile )
    tables = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name", conn)
    query = 'select * from %s'
    tables = { tname: pd.read_sql(query % tname, conn) for tname in tables.name.values }
    return tables

class cMonkey2:
    dbfile = None
    tables = None
    iteration = 2001
    k_clust = None

    def __init__( self, dbfile ):
        self.dbfile = dbfile
        self.tables = read_all_tables( dbfile )
        self.iteration = max(self.tables['motif_infos'].iteration)
        self.k_clust = self.tables['run_infos'].num_clusters[0] ##max(self.tables['row_members'].cluster)

    def get_rows( self, k ):
        t1 = self.tables['row_members']
        t1 = t1[ (t1.iteration == self.iteration) & (t1.cluster == k) ]
        t2 = self.tables['row_names']
        t2 = pd.merge( t1, t2, on='order_num' )
        return t2.name.values

    def get_cols( self, k ):
        t1 = self.tables['column_members']
        t1 = t1[ (t1.iteration == self.iteration) & (t1.cluster == k) ]
        t2 = self.tables['column_names']
        t2 = pd.merge( t1, t2, on='order_num' )
        return t2.name.values

    def get_cluster_info( self, k ):
        t1 = self.tables['cluster_stats']
        t1 = t1[ t1.cluster == k ]
        t1 = t1.drop( ['iteration', 'cluster'], 1 )

        t2 = self.tables['motif_infos']
        t2 = t2[ t2.cluster == k ]
        t2 = t2.drop( ['iteration', 'cluster'], 1 )
        return {'cluster':t1, 'motif':t2}

    def clusters_w_genes( self, genes ):
        t1 = self.tables['row_members']
        t1 = t1[ (t1.iteration == self.iteration) ]
        t2 = self.tables['row_names']
        t2 = t2[ np.in1d(t2.name, genes) ]
        t2 = pd.merge( t1, t2, on='order_num' )
        t2 = t2.drop( ['iteration', 'order_num'], 1 )
        return t2

    def clusters_w_conds( self, conds ):
        t1 = self.tables['column_members']
        t1 = t1[ (t1.iteration == self.iteration) ]
        t2 = self.tables['column_names']
        t2 = t2[ np.in1d(t2.name, conds) ]
        t2 = pd.merge( t1, t2, on='order_num' )
        t2 = t2.drop( ['iteration', 'order_num'], 1 )
        return t2

    def get_motif_pssm(self, cluster_num, motif_num):
        """export the specified motif to a pandas dataframe
        Parameters:
        - cluster_num: bicluster number
        - motif_num: motif number
        """
        #conn = sqlite3.connect(self.dbfile)
        #cursor = conn.cursor()
        #cursor.execute('select max(iteration) from motif_infos')
        #iteration = cursor.fetchone()[0]

        #query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
        #params = [self.iteration, cluster_num, motif_num]
        #cursor.execute(query, params)
        #rowid = cursor.fetchone()[0]
        motif_infos = self.tables['motif_infos']
        rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
                            (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1

        #query = 'select a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?'
        #params = [self.iteration, rowid]
        #pssm = pd.read_sql( query, conn, params=params )
        motif_pssm_rows = self.tables['motif_pssm_rows']
        pssm = motif_pssm_rows[(motif_pssm_rows.iteration==self.iteration) & (motif_pssm_rows.motif_info_id==rowid)]
        pssm.drop( ['motif_info_id', 'iteration', 'row'], 1, inplace=True )
        return pssm

    def get_biop_motif(self, cluster_num, motif_num, option='sites'):
        ##import egrin2.export_motifs as em
        """export the specified motif to a biopython motif object
        Parameters:
        - cluster_num: bicluster number
        - motif_num: motif number
        - option of how to translate - sites: jaspar 'sites' file; pfm: jaspar 'pfm' file
        """
        #conn = sqlite3.connect(self.dbfile)
        #cursor = conn.cursor()
        #cursor.execute('select max(iteration) from motif_infos')
        #iteration = cursor.fetchone()[0]

        #query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
        #params = [self.iteration, cluster_num, motif_num]
        #cursor.execute(query, params)
        #rowid = cursor.fetchone()[0]
        motif_infos = self.tables['motif_infos']
        rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
                            (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        #mot_info = pd.read_sql('select * from motif_infos where rowid=?', conn, params=[rowid])
        #mot_sites = pd.read_sql('select * from meme_motif_sites where motif_info_id=?', conn, params=[rowid])
        mot_sites = self.tables['meme_motif_sites'][self.tables['meme_motif_sites'].motif_info_id == rowid]
            
        output = cStringIO.StringIO()
        ## ONE WAY TO TRY -- but Bio.motifs cant parse the incomplete MEME file
        ##output.write(em.MEME_FILE_HEADER % (0.25, 0.25, 0.25, 0.25))
        ##em.write_pssm(output, cursor, os.path.dirname(self.dbfile), cluster_num, rowid,
        ##              motif_num, mot_info['evalue'][0], 10)
        ##output.seek(0)
        ##mot = motifs.read( output, 'meme' )
            
        ## Second way - create a jaspar 'pfm' file from the pssm
        if option == 'pfm':
            #query = 'select a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?'
            #params = [self.iteration, rowid]
            #pssm = pd.read_sql( query, conn, params=params )

            motif_pssm_rows = self.tables['motif_pssm_rows']
            pssm = motif_pssm_rows[(motif_pssm_rows.iteration==self.iteration) & (motif_pssm_rows.motif_info_id==rowid)]
            pssm = pssm.drop( ['motif_info_id', 'iteration', 'row'], 1 )

            counts = np.round( pssm * mot_sites.shape[0] ).transpose()
            counts.to_string(output, header=False, index=False )
            output.seek(0)
            mot = motifs.read( output, 'pfm' )

            ## Third way - create a jaspar 'sites' file
        elif option == 'sites':
            seqs = {}
            for i in mot_sites.index.values:
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
        ## Note Bio.motifs.weblogo() uses the weblogo server (slow? requires connection.)
        #kwargs = dict(color_scheme='classic')
        #mot.weblogo('file.png', color_scheme='color_classic') ## note, can use format='PDF'
        #img = mpimg.imread('file.png')
        #imgplot = plt.imshow( img )
        #plt.show()
        return mot

    ## This uses weblogolib package to create files directly (installed as weblogo via pip)
    ## https://code.google.com/p/weblogo/
    def plot_motif( self, cluster_num, motif_num ):
        #conn = sqlite3.connect(self.dbfile)
        #cursor = conn.cursor()
        #cursor.execute('select max(iteration) from motif_infos')
        #iteration = cursor.fetchone()[0]

        #query = 'select rowid from motif_infos where iteration=? and cluster=? and motif_num=?'
        #params = [self.iteration, cluster_num, motif_num]
        #cursor.execute(query, params)
        #rowid = cursor.fetchone()[0]
        #mot_info = pd.read_sql('select * from motif_infos where rowid=?', conn, params=[rowid])
        #mot_sites = pd.read_sql('select * from meme_motif_sites where motif_info_id=?', conn, params=[rowid])

        motif_infos = self.tables['motif_infos']
        rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
                            (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        mot_sites = self.tables['meme_motif_sites'][self.tables['meme_motif_sites'].motif_info_id == rowid]

        ldata = wl.LogoData.from_seqs(wl.SeqList(mot_sites.seq.values.tolist(), wl.unambiguous_dna_alphabet))
        options = wl.LogoOptions()
        options.fineprint = os.path.dirname(self.dbfile) + ' %03d %03d' % ( cluster_num, motif_num )
        format = wl.LogoFormat(ldata, options) 
        format.color_scheme = wl.classic
        format.resolution = 150
        tmp = wl.png_formatter( ldata, format )
        
        output = cStringIO.StringIO(tmp)
        img = mpimg.imread(output)
        imgplot = plt.imshow( img )
        plt.show()
        return img


##from egrin2.cmonkeyobj import cMonkey2 as cm2
##b = cm2('eco-out-001/cmonkey_run.db')
