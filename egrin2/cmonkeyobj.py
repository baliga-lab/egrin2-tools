import numpy as np
import os
import pandas as pd
import cStringIO
import sqlite3
import bz2
import cPickle as pickle
import ConfigParser

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
    organism = None
    species = None
    ratios = None
    config = None

    def __init__( self, dbfile ):
        self.dbfile = dbfile
        self.tables = read_all_tables( dbfile )
        self.iteration = max(self.tables['motif_infos'].iteration)
        self.k_clust = self.tables['run_infos'].num_clusters[0] ##max(self.tables['row_members'].cluster)
        self.organism = self.tables['run_infos'].organism[0]
        self.species = self.tables['run_infos'].species[0]
        self.config = self.load_config()

    def get_feature_names( self ):
        feature_names_file = './cache/' + self.species + '_feature_names'
        feature_names = pd.read_table( feature_names_file, sep='\t', header=None, skiprows=4 )
        feature_names.columns = ['id','names','type']
        #feature_names = feature_names.set_index( 'names' )
        return feature_names

    def get_features( self ):
        features_file = './cache/' + self.species + '_features'
        features = pd.read_table( features_file, sep='\t', header=None, skiprows=16 )
        cols = features.columns.values; cols[0] = 'id'; features.columns = cols
        #features = features.set_index( 'od' )
        return features

    def load_ratios( self, ratios_file=None ):
        import gzip
        if ratios_file is None:
            ratios_file = os.path.dirname(self.dbfile) + '/ratios.tsv.gz'
        self.ratios = pd.read_table( gzip.GzipFile( ratios_file ), sep='\t' )
        return self.ratios

    def load_config( self, config_file=None ):
        """then can do e.g., b.config.getfloat('Rows', 'scaling_constant')
           or simply, dict(b.config.items('Rows'))"""
        if config_file is None:
            config_file = os.path.dirname(self.dbfile) + '/final.ini'
        config_parser = ConfigParser.ConfigParser()
        config_parser.read( config_file )
        self.config = config_parser
        return self.config

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

    def get_ratios( self, k=None, rows=None, cols=None, included=True ):
        """Extract submatrix of ratios for cluster or rows/cols. 
        If ~included, extract submatrix of ratios for conditions NOT in cluster."""
        if self.ratios is None:
            ratios = self.load_ratios()
        if k is not None:
            if rows is None:
                rows = self.get_rows( k )
            if cols is None:
                cols = self.get_cols( k )
        if not included:
            cols = ratios.columns.values[ np.in1d( ratios.columns.values, cols, invert=True ) ]
        rats = self.ratios.ix[ rows, cols ]
        return rats

    def plot_ratios( self, k=None, rows=None, cols=None, included=True, kind='line' ):
        ## see http://pandas.pydata.org/pandas-docs/version/0.15.0/visualization.html -- cool!
        ## can use kind = 'box' too!
        rats = self.get_ratios( k, rows, cols, included )
        rats = rats.transpose()

        if kind == 'box': ## sort by mean of columns
            means = rats.mean(1)
            tmp = pd.concat( [rats, means], 1 )
            cols = tmp.columns.values; cols[-1] = 'MEANS'; tmp.columns = cols
            tmp = tmp.sort( ['MEANS'] )
            tmp = tmp.drop( 'MEANS', 1 )
            rats = tmp.transpose()
            rats.plot(kind=kind, use_index=False, title='Cluster %d'%(k), legend=False, sym='.')
        else:
            rats.plot(kind=kind, use_index=False, title='Cluster %d'%(k), legend=False)
        ## use plt.close() to close the window

    def get_cluster_info( self, k ):
        t1 = self.tables['cluster_stats']
        t1 = t1[ t1.cluster == k ]
        #t1 = t1.drop( ['iteration', 'cluster'], 1 )

        t2 = self.tables['motif_infos']
        t2 = t2[ t2.cluster == k ]
        #t2 = t2.drop( ['iteration', 'cluster'], 1 )

        ## Extract it.
        out = {'residual':t1.residual.values[0],
               'nrows':t1.num_rows.values[0],
               'ncols':t1.num_cols.values[0],
               'e_values:':t2.evalue.values}

        ## Also get p-clust
        pclusts = np.array([self.get_motif_pclust(k,i) for i in range(1,t2.shape[0]+1)])
        out['pclusts'] = pclusts
        
        return out

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

    def __get_motif_id(self, cluster_num, motif_num):
        motif_infos = self.tables['motif_infos']
        rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
                            (motif_infos.cluster==cluster_num) & 
                            (motif_infos.motif_num==motif_num)].index.values[0]+1
        return rowid

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
        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)

        #query = 'select a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?'
        #params = [self.iteration, rowid]
        #pssm = pd.read_sql( query, conn, params=params )
        motif_pssm_rows = self.tables['motif_pssm_rows']
        pssm = motif_pssm_rows[(motif_pssm_rows.iteration==self.iteration) & (motif_pssm_rows.motif_info_id==rowid)]
        pssm.drop( ['motif_info_id', 'iteration', 'row'], 1, inplace=True )
        return pssm

    def get_motif_sites(self, cluster_num, motif_num):
        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)

        sites = self.tables['meme_motif_sites']
        sites = sites[ sites.motif_info_id == rowid ]
        sites = sites.drop( ['motif_info_id'], 1 )

        feature_names = self.get_feature_names()
        tmp = pd.merge( sites, feature_names, left_on='seq_name', right_on='id' )
        tmp = tmp[ np.in1d( tmp.names.values, self.tables['row_names'].name.values ) ]
        tmp = tmp.drop( ['seq_name', 'type'], 1 )
        tmp = tmp.drop_duplicates()

        return tmp ## need to update genes based on synonyms

    def get_motif_pclust(self, cluster_num, motif_num):
        rowid = self.__get_motif_id(cluster_num, motif_num)
        sites = self.tables['meme_motif_sites']
        sites = sites[ sites.motif_info_id == rowid ]
        #sites = sites.drop( ['motif_info_id'], 1 )
        return np.mean( np.log10(sites.pvalue.values) )

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
        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)
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

        #motif_infos = self.tables['motif_infos']
        #rowid = motif_infos[(motif_infos.iteration==self.iteration) & 
        #                    (motif_infos.cluster==cluster_num) & (motif_infos.motif_num==motif_num)].index.values[0]+1
        rowid = self.__get_motif_id(cluster_num, motif_num)
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
        plt.axis('off')
        imgplot = plt.imshow( img )
        #plt.show()
        return plt


##from egrin2.cmonkeyobj import cMonkey2 as cm2
##b = cm2('eco-out-001/cmonkey_run.db')
##pd.Series([b.get_cluster_info(k)['residual'] for k in range(1,b.k_clust)]).plot(kind='hist',bins=20)
##pd.DataFrame([b.get_cluster_info(k)['pclusts'] for k in range(1,b.k_clust)]).plot(kind='hist',bins=20,stacked=True)
