#!/usr/bin/python

import bz2
import glob
import optparse
import numpy as np
#import scipy.sparse as sparse
import os
import pandas as pd
import igraph as ig
import itertools
import cPickle as pickle   ## make shelve much faster!
import shelve ## shove    ## see http://www.evilchuck.com/2008/02/tell-python-to-shove-it.html about shove

def system( cmd ):
    print cmd
    tmp = os.popen( cmd ).read()
    return tmp

op = optparse.OptionParser()
op.add_option('-i', '--input_dir', default='tomtom_out', help="The location of tomtom results, bzip'd")
op.add_option('-opt', '--option', default=2, help="Filtering option, 1 or 2; see comments")
opt, args = op.parse_args()
if not opt.input_dir:
    op.error('need --input_dir option.  Use -h for help.')

## if option == 1, then we SUM all the log-pvals for multiple occurrences of the same motif pair, 
## then filter the sums to only those that are < -10 (same as done for original Halo ensemble)
## if option == 2, then we take the LOWEST log-pval over multiple occurrences of the same motif 
## pair, with no additional filtering
option = 2

input_dir = 'tomtom_out'
input_dir = opt.input_dir
files = glob.glob( input_dir + "/*tomtom.tsv.bz2" ) # folder with the tomtom files bzip'd
##dfs = {}
## can pd.concat work on shelved dataframes? YES. Note protocol=2 is faster and smaller.
dfs = shelve.open(input_dir+'/shelf.db', protocol=2, writeback=False)
##dfs = shove.Shove('sqlite:///'+input_dir+'/shove.db', compress=True) ## note this requires SQLAlchemy installed
##dfs = shove.Shove('file://./'+input_dir+'/shove.db', compress=True) 

if len(dfs) != len(files): ## if using a shelf, once this is done once, you don't have to do it again.
    for f in files:    
        print(f)
        if f in dfs.keys():
            continue
        try:
            df = pd.read_table( bz2.BZ2File(f), sep='\t' )
            print df.shape
            df.drop('Query consensus', 1, inplace=True)
            df.drop('Target consensus', 1, inplace=True)
            df = df[ df['p-value'] <= 0.01 ]
            df = df[ df['#Query ID'] != df['Target ID'] ]
            print df.shape
            dfs[ f ] = df
        except:
            continue

        #if f == files[20]:
        #   break

dfs2 = pd.concat( dfs.values(), axis=0 )  ## works with 'shelve's but not (compressed) 'shove's
dfs.close()
print dfs2.shape

if option == 2:
    dfs2 = dfs2.sort( 'p-value' ) ## sort so lower p-values come first (these are kept by drop_duplicates)
    dfs2 = dfs2.drop_duplicates( ['#Query ID', 'Target ID'] ) ## no, we sum up the duplicate weights below
    print dfs2.shape

## incase we fail on steps below...
#dfs2.to_csv( bz2.BZ2File('motifs_tomtom.tsv.bz2', 'w'), sep='\t', index=False, header=True )
#dfs2 = pd.read_table( bz2.BZ2File('motifs_tomtom.tsv.bz2', 'r') )

all_motifs = np.sort( np.unique( np.concatenate( (np.unique( dfs2['#Query ID'].values ), 
                                                  np.unique( dfs2['Target ID'].values )), 0 ) ) )
all_motifs_lookup = dict( itertools.izip( all_motifs.tolist(), range(len(all_motifs)) ) )
inds1 = [ all_motifs_lookup[ i ] for i in dfs2['#Query ID'].values ]
inds2 = [ all_motifs_lookup[ i ] for i in dfs2['Target ID'].values ]

#X = sparse.lil_matrix( (len(all_motifs), len(all_motifs)) )
#X[ inds1, inds2 ] = -np.log10( dfs['p-value'].values )

#gr = pd.DataFrame( {'query':inds1, 'target':inds2, 'weight':np.round_(-np.log10(dfs2['p-value'].values+1e-99), 4)} )
#gr.to_csv( bz2.BZ2File('motifs_graph.tsv.bz2', 'w'), sep=' ', index=False, header=False )
#gr2 = ig.Graph.Read_Ncol( bz2.BZ2File('motifs_graph.tsv.bz2'), names=False, weights=True, directed=False )
gr2 = ig.Graph( zip(inds1, inds2), edge_attrs={'weight':np.round_(-np.log10(dfs2['p-value'].values+1e-99), 4)} )
print gr2.ecount(), gr2.vcount()

## cool! see http://igraph.org/python/doc/igraph.GraphBase-class.html#simplify
## add up the weights for duplicated edges into a single edge weight
if option == 1:
    gr2a = gr2.simplify( multiple=True, loops=False, combine_edges=sum ) 
    print gr2a.ecount(), gr2a.vcount()
    ## see http://igraph.org/python/doc/tutorial/tutorial.html#selecting-vertices-and-edges
    gr2b = gr2a.es.select( weight_gt=10 )  ## returns an EdgeList
    gr2b = gr2b.subgraph()   ## convert to a graph
    print gr2b.ecount(), gr2b.vcount()
elif option == 2:
    gr2a = gr2.simplify( multiple=True, loops=False, combine_edges=max )
    print gr2a.ecount(), gr2a.vcount()
    gr2b = gr2a

gr2b.write_ncol( "mot_metaclustering.txt" )
##gr2b = ig.Graph.Read_Ncol( 'mot_metaclustering.txt', names=True, weights=True, directed=False )    

## now run mcl, latest version from http://www.micans.org/mcl/src/mcl-latest.tar.gz
system( "../mcl mot_metaclustering.txt --abc -I 4.5 -v all -te 3 -S 200000" )

fo = open( 'out.mot_metaclustering.txt.I45', 'r' )
lines = fo.readlines()
fo.close()
lines = [np.array(line.split(), np.int) for line in lines]
clusters = [all_motifs[i] for i in lines]

lines2 = ['\t'.join(i) for i in clusters]
fo = open( 'out.mot_metaclustering.txt.I45.txt', 'w' )
fo.write( '\n'.join(lines2) )
fo.close()

clust_lens = np.array([len(i) for i in clusters])
print 'Clusters with >= 10 motifs:', np.sum( clust_lens >= 10 )
print 'Total number of motifs:', gr2a.vcount()
print 'Number of motifs in >= 10-size clusters:', \
    np.sum( np.array( [ len(clusters[i]) for i in np.where( clust_lens >= 10 )[0] ] ) )
print 'Fraction of motifs in >= 10-size clusters:', \
    float( np.sum( np.array( [ len(clusters[i]) for i in np.where( clust_lens >= 10 )[0] ] ) ) ) / \
    float( gr2a.vcount() )

## Get info on alignments for each motif cluster
dfs2 = dfs2.set_index('#Query ID', drop=False)
clust_dfs = {}
for i in xrange( len(clusters) ):
    clust = clusters[i]
    print i, len(clust)
    if i in clust_dfs.keys() or len(clust) < 10 or i > 500:
        continue
    df = dfs2.ix[ clust ] ##[ np.in1d( dfs2['#Query ID'].values, clust ) ]
    df = df.loc[ np.in1d( df['Target ID'].values, clust ) ]
    print df.shape
    clust_dfs[i] = df

## Create the "aligned pssms" and "combined pssm" for each motif cluster
## Not ready yet.
# for i in xrange( len(clusters) ):
#     df = clust_dfs[i]
#     df = df.sort( 'p-value' )

#     all_pssms = {}
#     for mot in clusters[i]:
#         tmp = np.array( mot.split('-')[2].split('_'), np.int )
#         dbfile = 'eco-out-%03d/cmonkey_run.db' % (tmp[0])
#         ## if use get_biop_motif(), then can do e.g., utils.plot_pssm(all_pssms['eco-out-009_312_02'])
#         ##all_pssms[mot] = utils.get_biop_motif( dbfile, tmp[1], tmp[2], 'pfm' ) ## sites opt. doesnt work well w 'N's
#         ## or can do utils.plot_motif(dbfile,tmp[1],tmp[2])
#         all_pssms[mot] = utils.get_motif_pssm( dbfile, tmp[1], tmp[2] )

#     widths = np.array( [ pssm.shape[0] for pssm in all_pssms.values() ] )
#     max_width = np.max( widths )
#     tmp = np.zeros(max_width)
#     empty_pssm = pd.DataFrame( { 'a':tmp, 'c':tmp, 'g':tmp, 't':tmp } )
    
#     aligned_pssms = {}

#     query = df.ix[0, '#Query ID']
#     qpssm = all_pssms[query]
#     tmp = empty_pssm.copy()

#     empty = pd.Series( { 'a':0., 'c':0., 'g':0., 't':0. } )
#     for j in xrange(df.shape[0]):
        
#         tmp = np.array( query.split('-')[2].split('_'), np.int )
#         dbfile = 'eco-out-%03d/cmonkey_run.db' % (tmp[0])
#         qpssm = get_motif_pssm( dbfile, tmp[1], tmp[2] )
#         if aligned_pssm is None: ## make empty data frame
#             aligned_pssm = qpssm.copy() * 0
#             while aligned_pssm.shape[0] < 64:
#                 aligned_pssm = pd.concat( [ aligned_pssm, aligned_pssm ] )
#             aligned_pssm = aligned_pssm[ 0:64 ]
#             aligned_pssm = aligned_pssm.set_index( np.arange(64) )

#         if not query in done:
#             aligned_pssm

#         target = df.ix[j, 'Target ID']
#         row = 

#     for mot in clusters[i]:
#         tmp = np.array( mot.split('-')[2].split('_'), np.int )
#         dbfile = 'eco-out-%03d/cmonkey_run.db' % (tmp[0])
#         pssm = get_motif_pssm( dbfile, tmp[1], tmp[2] )
        
