#!/usr/bin/python

import os
import bz2
import glob
import itertools
#import cPickle as pickle   ## make shelve much faster!
#import shelve ## shove    ## see http://www.evilchuck.com/2008/02/tell-python-to-shove-it.html about shove
import optparse
import numpy as np
import numpy.core.defchararray as npstr
#import scipy.sparse as sparse
import pandas as pd
import igraph as ig

def system( cmd ):
    print cmd
    tmp = os.popen( cmd ).read()
    return tmp

op = optparse.OptionParser()
op.add_option('-i', '--input_dir', default='tomtom_out', help="The location of tomtom results, bzip'd")
op.add_option('-f', '--filter_motifs', default=False, action='store_true', 
              help="Pre-filter motifs to include; see comments")
op.add_option('--plot_motifs', default=False, action='store_true', 
              help="Plot motif clusters in separate PDFs")
op.add_option('-o', '--option', default=1, help="Filtering option, 1 or 2; see comments")
opt, args = op.parse_args()
if not opt.input_dir:
    op.error('need --input_dir option.  Use -h for help.')

## if option == 1, then we SUM all the log-pvals for multiple occurrences of the same motif pair, 
## then filter the sums to only those that are < -10 (same as done for original Halo ensemble)
## if option == 2, then we take the LOWEST log-pval over multiple occurrences of the same motif 
## pair, with no additional filtering. Note, option (2) was used for Eco and (1) for Halo in MSB EGRIN2 paper.
option = 1
if opt.option:
    option = opt.option

print 'OPTION:', option

## pre-filter motifs, not implemented yet. (1) remove motifs that are in coding regions (from fimo table);
## (2) filter by motif E-value (3) filter by bicluster residual?
pre_filter = False
if opt.filter_motifs:
    pre_filter = opt.filter_motifs

coding_fracs = total_frac = None

if pre_filter:
    try: ## necessary only because egrin2-tools has hyphen and can't have hyphens in python module paths...
        os.symlink('egrin2-tools/postproc/coding_fracs.py', 'coding_fracs.py')
    except:
        None
    import coding_fracs as cf
    total_frac = cf.get_total_coding_rgn('cache/Mycobacterium_tuberculosis_H37Rv_features')

    cf_files = np.sort( np.array( glob.glob(os.path.join('*/coding_fracs.tsv.bz2')) ) )
    coding_fracs = [] ##{}
    for f in cf_files:
        print f
        cff = pd.read_table( bz2.BZ2File(f), sep='\t' )
        cm_run = os.path.dirname(f) ##.split('-')[2]
        cff['cm_run'] = cm_run
        if cff.shape[0] > 1:
            coding_fracs.append( cff ) ##[f] = cff 
    coding_fracs = pd.concat( coding_fracs, keys=None, ignore_index=True )
    ## this has a hack - for some reason cluster_id in coding_fracs is %04d, trim first zero to make it %03d ...
    splitted = npstr.split( coding_fracs.motif.values.astype(str), '_' )
    clust_id = np.char.mod('_%03d_', np.array( [int(i[0]) for i in splitted] ) ) ## see https://stackoverflow.com/a/28286749
    mot_id = np.char.mod('%02d', np.array( [int(i[1]) for i in splitted] ) )
    mot_id = npstr.add( clust_id, mot_id )
    coding_fracs[ 'motif_id' ] = npstr.add( coding_fracs.cm_run.values.astype(str), mot_id )

input_dir = 'tomtom_out'
input_dir = opt.input_dir
files = np.sort( np.array( glob.glob( input_dir + "/*tomtom.tsv.bz2" ) ) ) # folder with the tomtom files bzip'd
dfs = {}
## can pd.concat work on shelved dataframes? YES. Note protocol=2 is faster and smaller.
#dfs = shelve.open('tomtom_shelf.db', protocol=2, writeback=False)
##dfs = shove.Shove('sqlite:///'+input_dir+'/shove.db', compress=True) ## note this requires SQLAlchemy installed
##dfs = shove.Shove('file://./'+input_dir+'/shove.db', compress=True) 

if len(dfs) != len(files): ## if using a shelf, once this is done once, you don't have to do it again.
    for f in files:    
        gene = os.path.basename(f).split('.')[0]
        print f, gene
        if gene in dfs.keys():
            continue
        try:
            df = pd.read_table( bz2.BZ2File(f), sep='\t' )
            print df.shape
            if df.shape[0] <= 0:
                continue
            ##df[ 'gene' ] = gene
            ##df = pd.concat( [ df, pd.Series(np.repeat(gene,df.shape[0])) ], 1 )
            #cols = df.columns.values; cols[-1] = 'gene'; df.columns = cols
            df = df.ix[ df['p-value'] <= 0.01 ] ## 0.1 ]
            print df.shape
            df = df.ix[ df['#Query ID'] != df['Target ID'] ]
            print df.shape
            df = df.ix[ df.Overlap >= 6 ] ## same setting as Halo run
            df = df.drop(['Query consensus', 'Target consensus'], axis=1)

            if pre_filter:
                ## add the coding fracs to the df:
                tmp = pd.merge(df, coding_fracs, how='left', left_on='#Query ID', right_on='motif_id')
                tmp = pd.merge(tmp, coding_fracs, how='left', left_on='Target ID', right_on='motif_id')
                tmp = tmp.drop( ['motif_x', 'cm_run_x', 'motif_id_x', 'motif_y', 'cm_run_y', 'motif_id_y'], axis=1 )
                ## drop the motifs with coding fracs greater than (expected value) + (obs. stddev)/2
                cutoff = total_frac ## + coding_fracs.coding_frac.mad() / 2
                tmp = tmp.ix[ np.logical_or( tmp.coding_frac_x.values < cutoff, tmp.coding_frac_y.values < cutoff ) ]
                print tmp.shape
                df = tmp; del tmp

            dfs[ gene ] = df

        except:
            continue

        #if f == files[20]:
        #   break

##dfs2 = pd.concat( dfs.values(), axis=0 )  ## works with 'shelve's but not (compressed) 'shove's
dfs2 = pd.concat( dfs, axis=0 )
#dfs.close()
print dfs2.shape

## incase we fail on steps below...
dfs2.to_csv( bz2.BZ2File('motifs_tomtom.tsv.bz2', 'w'), sep='\t', index=False, header=True )
#dfs2 = pd.read_table( bz2.BZ2File('motifs_tomtom.tsv.bz2', 'r') )

if option == 2:
    dfs2.sort( 'p-value', inplace=True ) ## sort so lower p-values come first (these are kept by drop_duplicates)
    dfs2.drop_duplicates( ['#Query ID', 'Target ID'], inplace=True ) ## no, we sum up the duplicate weights below
    print dfs2.shape

# all_motifs = np.sort( np.unique( np.concatenate( (np.unique( dfs2['#Query ID'].values ), 
#                                                   np.unique( dfs2['Target ID'].values )), 0 ) ) )
# all_motifs_lookup = dict( itertools.izip( all_motifs.tolist(), range(len(all_motifs)) ) )
# inds1 = [ all_motifs_lookup[ i ] for i in dfs2['#Query ID'].values ]
# inds2 = [ all_motifs_lookup[ i ] for i in dfs2['Target ID'].values ]

#X = sparse.lil_matrix( (len(all_motifs), len(all_motifs)) )
#X[ inds1, inds2 ] = -np.log10( dfs['p-value'].values )

gr = pd.DataFrame( {'query':dfs2['#Query ID'].values, 'target':dfs2['Target ID'].values, 
                    'weight':np.round_(-np.log10(dfs2['p-value'].values+1e-99), 4)} )
#gr = pd.DataFrame( {'query':inds1, 'target':inds2, 'weight':np.round_(-np.log10(dfs2['p-value'].values+1e-99), 4)} )
gr.to_csv( 'motifs_graph.tsv', sep=' ', index=False, header=False ) ## igraph cannot read bzipped files (streams)
del gr
gr2 = ig.Graph.Read_Ncol( 'motifs_graph.tsv', names=True, weights=True, directed=False )
ut.system( 'pbzip2 -fv9 motifs_graph.tsv &' )
#gr2 = ig.Graph( zip(inds1, inds2), edge_attrs={'weight':np.round_(-np.log10(dfs2['p-value'].values+1e-99), 4)} )
print gr2.ecount(), gr2.vcount()
#gr2.write_ncol( "mot_metaclustering.txt" )

## cool! see http://igraph.org/python/doc/igraph.GraphBase-class.html#simplify
## add up the weights for duplicated edges into a single edge weight
if option == 1:
    gr2a = gr2.simplify( multiple=True, loops=False, combine_edges=sum ) 
    print gr2a.ecount(), gr2a.vcount()
    ## see http://igraph.org/python/doc/tutorial/tutorial.html#selecting-vertices-and-edges
    gr2b = gr2a.es.select( weight_gt=10 )  ## used 10 for Halo; use less for fewer runs. returns an EdgeList
    gr2b = gr2b.subgraph()   ## convert to a graph
    print gr2b.ecount(), gr2b.vcount()
elif option == 2:
    gr2a = gr2.simplify( multiple=True, loops=False, combine_edges=max )
    print gr2a.ecount(), gr2a.vcount()
    gr2b = gr2a

del gr2
gr2b.write_ncol( "mot_metaclustering.txt", weights=None ) ## no weights used - same as Halo analysis which was best!
##gr2b = ig.Graph.Read_Ncol( 'mot_metaclustering.txt', names=True, weights=True, directed=False )    

## now run mcl, latest version from http://www.micans.org/mcl/src/mcl-latest.tar.gz
system( "./progs/mcl mot_metaclustering.txt --abc -I 4.5 -v all -te 3 -S 200000" )

fo = open( 'out.mot_metaclustering.txt.I45', 'r' )
lines = fo.readlines()
fo.close()
##lines = [np.array(line.split(), np.int) for line in lines]
lines = [np.array(line.split()) for line in lines]
clusters = lines   ## file contains actual motif ids rather than numbers
del lines

# clusters = [all_motifs[i] for i in lines]
# lines2 = ['\t'.join(i) for i in clusters]
# fo = open( 'out.mot_metaclustering.txt.I45.txt', 'w' )
# fo.write( '\n'.join(lines2) )
# fo.close()

clust_lens = np.array([len(i) for i in clusters])
print 'Clusters with >= 10 motifs:', np.sum( clust_lens >= 10 )
print 'Total number of motifs:', gr2a.vcount()
print 'Number of motifs in >= 10-size clusters:', \
    np.sum( np.array( [ len(clusters[i]) for i in np.where( clust_lens >= 10 )[0] ] ) )
print 'Fraction of motifs in >= 10-size clusters:', \
    float( np.sum( np.array( [ len(clusters[i]) for i in np.where( clust_lens >= 10 )[0] ] ) ) ) / \
    float( gr2a.vcount() )
del gr2a

## Get info on alignments for each motif cluster
dfs3 = dfs2.set_index('#Query ID', drop=False)
clust_dfs = {}
for i in xrange( len(clusters) ):
    clust = clusters[i]
    print i, len(clust)
    if i in clust_dfs.keys() or len(clust) < 10 or i > 500:
        continue
    df = dfs3.ix[ clust ] ##[ np.in1d( dfs3['#Query ID'].values, clust ) ]
    df = df.iloc[ np.in1d( df['Target ID'].values, clust ) ]
    df = df.sort( ['p-value'] )
    df = df.ix[ ~df.duplicated(['#Query ID','Target ID']) ] ## remove dupes
    df = df.reset_index( drop=True )
    df[ 'motif_clust' ] = i
    print df.shape
    clust_dfs[i] = df

del dfs3
clust_dfs = pd.concat( clust_dfs, axis=0 )
print clust_dfs.shape
## get coding fracs per motif cluster via:
## clust_dfs.groupby('motif_clust').mean().coding_frac_x

clust_dfs.to_csv( bz2.BZ2File('motif_clusts.tsv.bz2', 'w'), sep='\t', index=False, header=True )
#clust_dfs = pd.read_table( bz2.BZ2File('motif_clusts.tsv.bz2', 'r') )

## Create the "aligned pssms" and(or?) "combined pssm" for each motif cluster
## Not ready yet.

if False:
    try: ## necessary to import cmonkey.meme
        os.symlink('cmonkey2/cmonkey', 'cmonkey')
    except:
        None

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio import SeqIO
    #from Bio import motifs
    #import StringIO
    import weblogolib as wl
    import cmonkey.meme as meme ## use Wei-Ju's meme parser - much better than BioPython!
    import utils as ut

    try:
        os.mkdir("motif_clusters") ## location of output files
    except:
        None

    def writeFasta( seqs, fname ):
        tmp = [ SeqRecord(Seq(seqs[i], IUPAC.ExtendedIUPACDNA()), id=i, name=i, description=i) for i in seqs.keys() ]
        handle = open(fname, 'w')
        SeqIO.write(tmp, handle, "fasta")
        handle.close()

    def plot_motif_from_sites( sites, img_format='png', smallText=None ):
        ldata = wl.LogoData.from_seqs(wl.SeqList(sites, wl.unambiguous_dna_alphabet))
        options = wl.LogoOptions()
        if smallText is not None:
            options.fineprint = smallText ##os.path.dirname(self.dbfile) + ' ' + self.organism
        format = wl.LogoFormat(ldata, options) 
        format.color_scheme = wl.classic
        format.resolution = 150
        if img_format == 'png':
            tmp = wl.png_formatter( ldata, format )
            output = cStringIO.StringIO(tmp)
            img = mpimg.imread(output)
            plt.axis('off')
            imgplot = plt.imshow( img )
            return plt
        elif img_format == 'svg':
            tmp = wl.svg_formatter( ldata, format )
            return tmp
        elif img_format == 'pdf':
            tmp = wl.pdf_formatter( ldata, format )
            return tmp

    def plot_multiple_motifs( cluster_ind, max_plot=16, outPdf=None ):
        from matplotlib import pyplot as plt
        pdf = None
        if outPdf is not None:
            from matplotlib.backends.backend_pdf import PdfPages
            pdf = PdfPages(outPdf)

        clust = clusters[ cluster_ind ]
        out = []
        plotted = 0
        for c in np.sort(clust):
            print c
            tmp = c.split('_')
            cm_ind = tmp[0] + '/cmonkey_run.db'
            cl = int(tmp[1])
            mot = int(tmp[2])
            cm = cms[ cm_ind ]
            plt.subplot(4, 4, plotted+1)
            cm.plot_motif( cl, mot )
            plotted += 1
            if plotted >= max_plot:
                if outPdf is not None:
                    pdf.savefig()
                break
        if outPdf is not None:
            pdf.close()

    from cmonkeyobj import cMonkey2 as cm2
    ## pre-read in all cmonkey database objects for speed
    cmdbs = sorted( glob.glob( 'mtu-out-*/cmonkey_run.db' ) )
    cms = { db:cm2(db) for db in cmdbs }
    
    for iii in xrange( 0,len(clusters) ):
        ##df = clust_dfs.ix[iii]  # dont need for this?
        if len(clusters[iii]) < 10:
            break
        memeOut = None
        if not os.path.exists( 'motif_clusters/%04d_memeOut.txt'%(iii) ):
            all_pssms = {}
            all_sites = {} ## maybe use the actual seqs and run meme on them?
            for mot in clusters[iii]:
                print iii, mot, len(clusters[iii])
                tmp = np.array( mot.split('-')[2].split('_'), np.int )
                dbfile = 'mtu-out-%03d/cmonkey_run.db' % (tmp[0])
                cm = cms[ dbfile ] ##cm2( dbfile )
                ## if use get_biop_motif(), then can do e.g., utils.plot_pssm(all_pssms['eco-out-009_312_02'])
                ##all_pssms[mot] = utils.get_biop_motif( dbfile, tmp[1], tmp[2], 'pfm' ) ## doesnt work well w 'N's
                ## or can do cm.plot_motif(tmp[1],tmp[2])
                all_pssms[mot] = cm.get_motif_pssm( tmp[1], tmp[2] )
                sites = cm.get_motif_sites( tmp[1], tmp[2] )
                all_sites[mot] = sites ##[['flank_left','seq','flank_right']]

            sites = pd.concat( all_sites, axis=0 )
            #sites[ pd.isnull(sites) ] = ""
            sites['tmp'] = np.log10(sites.pvalue.values) / np.array([len(i) for i in sites.seq.values]).astype(float)
            sites.sort( ['tmp'], inplace=True )
            tmp = sites.groupby('names')
            firsts = tmp.first() ## get first row of each group
            sizes = tmp.size()
            firsts.fillna('', inplace=True)
            if len(clusters[iii]) > 100:
                firsts = firsts[ sizes > len(clusters[iii])/10. ]
            else:
                firsts = firsts[ sizes >= 2 ]
            maxlen = max( [len(i) for i in firsts.seq.values] )
            minlen = min( [len(i) for i in firsts.seq.values] )
            seqs = [ ( firsts.flank_left[i] + firsts.seq[i] + \
                                  firsts.flank_right[i] ).strip('X').replace('X','N') \
                                for i in range(firsts.shape[0]) ]
            print len(seqs)
            seqs = { firsts.index.values[i]: ('N'*20)+seqs[i].strip('X').replace('X','N')+('N'*20) \
                     for i in range(len(seqs)) }
            writeFasta( seqs, 'motif_clusters/%04d.fna'%(iii) )
            cmd = './progs/meme %s -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %d -evt 999999 ' \
                  '-minw %d -maxw %d -mod oops -nostatus -text -time 600' \
                  % ('motif_clusters/%04d.fna'%(iii), 1, minlen, maxlen)
            print cmd
            try:
                memeOut = os.popen( cmd ).read()
            except:
                print 'MEME DID NOT COMPLETE WITHIN 10 MINUTES'
                continue
        else:
            fo = open( 'motif_clusters/%04d_memeOut.txt'%(iii), 'r' )
            memeOut = fo.read()
            fo.close()
        #handle = StringIO.StringIO(memeOut)
        #record = motifs.parse(handle, "meme")
        #handle.close()
        try:
            record = meme.read_meme_output(memeOut, 1)[0]
        except:
            print 'COULD NOT PARSE MEME OUTPUT'
            continue
        
        sites = pd.DataFrame(record.sites)[5].values
        smallText = str( iii ) + ' ' + str(len(clusters[iii])) + ' ' + str(record.num_sites) + \
                    ' E=' + str(record.evalue) 
        pdf_data = plot_motif_from_sites( sites, 'pdf', smallText )
        ut.writeLines( pdf_data.split( '\n' ), 'motif_clusters/%04d.pdf'%(iii) )
        ut.writeLines( memeOut.split( '\n' ), 'motif_clusters/%04d_memeOut.txt'%(iii) )
        print iii, 'DONE'
        
    os.popen( 'pdftk motif_clusters/????.pdf cat output motif_clusters/ALL.pdf' ).read()
    os.popen( 'pdfnup --landscape --suffix nup --nup 6x6 motif_clusters/ALL.pdf' ).read()

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
        
