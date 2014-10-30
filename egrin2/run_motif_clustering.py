#!/usr/bin/python

import bz2
import glob
#import time
#import math
import numpy as np
#import scipy.sparse as sparse
#import os
import pandas as pd
import igraph as ig
import itertools

def system( cmd ):
    print cmd
    tmp = os.popen( cmd ).read()
    return tmp

def main():
    files = glob.glob("tomtom_out/*tomtom.tsv.bz2") # folder with the tomtom files bzip'd
    dfs = {}
    for f in files:
        print(f)
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
            break
        #if f == files[20]:
        #   break
        
    dfs = pd.concat( dfs.values(), axis=0 )

    dfs2 = dfs.sort( 'p-value' ) ## sort so lower p-values come first (these are kept by drop_duplicates)
    dfs2 = dfs2.drop_duplicates( ['#Query ID', 'Target ID'] )
    all_motifs = np.sort( np.unique( np.concatenate( (np.unique( dfs2['#Query ID'].values ), 
                                                      np.unique( dfs2['Target ID'].values )), 0 ) ) )
    all_motifs_lookup = dict( itertools.izip( all_motifs.tolist(), range(len(all_motifs)) ) )
    inds1 = [ all_motifs_lookup[ i ] for i in dfs['#Query ID'].values ]
    inds2 = [ all_motifs_lookup[ i ] for i in dfs['Target ID'].values ]

    #X = sparse.lil_matrix( (len(all_motifs), len(all_motifs)) )
    #X[ inds1, inds2 ] = -np.log10( dfs['p-value'].values )

    #gr = pd.DataFrame( {'query':inds1, 'target':inds2, 'weight':np.round_(-np.log10(dfs['p-value'].values+1e-99), 4)} )
    #gr.to_csv( bz2.BZ2File('motifs_graph.tsv.bz2', 'w'), sep=' ', index=False, header=False )
    #gr2 = ig.Graph.Read_Ncol( bz2.BZ2File('motifs_graph.tsv.bz2'), names=False, weights=True, directed=False )
    gr2 = ig.Graph( zip(inds1, inds2), edge_attrs={'weight':np.round_(-np.log10(dfs['p-value'].values+1e-99), 4)} )
    print gr2.ecount(), gr2.vcount()
    ## cool! see http://igraph.org/python/doc/igraph.GraphBase-class.html#simplify
    gr2a = gr2.simplify( loops=False, combine_edges=sum ) 
    print gr2a.ecount(), gr2a.vcount()

    ## see http://igraph.org/python/doc/tutorial/tutorial.html#selecting-vertices-and-edges
    gr2b = gr2a.es.select( weight_gt=10 )  ## returns an EdgeList
    gr2b = gr2b.subgraph()   ## convert to a graph
    print gr2b.ecount(), gr2b.vcount()
    
    gr2b.write_ncol( "mot_metaclustering.txt" )
    ##gr2b = ig.Graph.Read_Ncol( 'mot_metaclustering.txt', names=True, weights=True, directed=False )    

    ## now run mcl, latest version from http://www.micans.org/mcl/src/mcl-latest.tar.gz
    system( "../mcl mot_metaclustering.txt --abc -I 4.5 -v all -te 3 -S 200000" )

    fo = open( 'out.mot_metaclustering.txt.I45', 'r' )
    lines = fo.readlines()
    fo.close()
    lines = [np.array(line.split(), np.int) for line in lines]
    lines = [all_motifs[i] for i in lines]

    lines2 = ['\t'.join(i) for i in lines]
    fo = open( 'out.mot_metaclustering.txt.I45.txt', 'w' )
    fo.write( '\n'.join(lines2) )
    fo.close()

if __name__ == '__main__':
    #start_time = float(time.time())
    main()
    #print "--- %f seconds ---" % (float(time.time()) - start_time)
