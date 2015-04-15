import os
import glob
import cStringIO

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
#from Bio import motifs
import weblogolib as wl
import utils as ut

try: ## necessary to import cmonkey.meme
    os.symlink('cmonkey2/cmonkey', 'cmonkey')
except:
    None

import cmonkey.meme as meme ## use Wei-Ju's meme parser - much better than BioPython!

organism = 'mtu' ## 'eco' ##
print organism, param_I_str

try:
    os.mkdir("motif_clusters_%s"%(param_I_str)) ## location of output files
except:
    None

fo = open( 'out.mot_metaclustering.txt.I%s'%(param_I_str), 'r' )
lines = fo.readlines()
fo.close()
##lines = [np.array(line.split(), np.int) for line in lines]
lines = [np.array(line.split()) for line in lines]
clusters = lines   ## file contains actual motif ids rather than numbers
del lines

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

def get_motif_cluster_sites( cluster_ind, force=False ):
    iii = cluster_ind
    memeOut = None
    if force or not os.path.exists( 'motif_clusters%s/%04d_memeOut.txt'%(param_I_str,iii) ):
        all_pssms = {}
        all_sites = {} ## maybe use the actual seqs and run meme on them?
        for mot in clusters[iii]:
            print iii, mot, len(clusters[iii])
            tmp = np.array( mot.split('-')[2].split('_'), np.int )
            dbfile = '%s-out-%03d/cmonkey_run.db' % (organism,tmp[0])
            cm = cms[ dbfile ] ##cm2( dbfile )
            ## if use get_biop_motif(), then can do e.g., utils.plot_pssm(all_pssms['eco-out-009_312_02'])
            ##all_pssms[mot] = utils.get_biop_motif( dbfile, tmp[1], tmp[2], 'pfm' ) ## doesnt work well w 'N's
            ## or can do cm.plot_motif(tmp[1],tmp[2])
            all_pssms[mot] = cm.get_motif_pssm( tmp[1], tmp[2] )
            sites = cm.get_motif_sites( tmp[1], tmp[2] )
            all_sites[mot] = sites ##[['flank_left','seq','flank_right']]

        sites = pd.concat( all_sites, axis=0 )
        #sites[ pd.isnull(sites) ] = ""
        ##sites['tmp'] = np.log10(sites.pvalue.values) / np.array([len(i) for i in sites.seq.values]).astype(float)
        sites['tmp'] = np.log10(sites.pvalue.values).astype(float)
        sites.sort( ['tmp'], inplace=True )
        tmp = sites.groupby('names')
        firsts = tmp.first() ## get first row of each group
        sizes = tmp.size()
        firsts.fillna('', inplace=True)
        if len(clusters[iii]) > 100:
            firsts = firsts[ sizes > len(clusters[iii])/20. ]
        else:
            firsts = firsts[ sizes > 2 ]
        maxlen = max( [len(i) for i in firsts.seq.values] )
        minlen = min( [len(i) for i in firsts.seq.values] )
        seqs = [ ( firsts.flank_left[i] + firsts.seq[i] + \
                   firsts.flank_right[i] ).strip('X').replace('X','N') \
                 for i in range(firsts.shape[0]) ]
        print len(seqs)
        seqs = { firsts.index.values[i]: ('N'*20)+seqs[i].strip('X').replace('X','N')+('N'*20) \
                 for i in range(len(seqs)) }
        writeFasta( seqs, 'motif_clusters_%s/%04d.fna'%(param_I_str,iii) )
        cmd = './progs/meme %s -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %d -evt 999999 ' \
              '-minw %d -maxw %d -mod oops -nostatus -text -time 600' \
              % ('motif_clusters_%s/%04d.fna'%(param_I_str,iii), 1, minlen, maxlen)
        print cmd
        try:
            memeOut = os.popen( cmd ).read()
        except:
            print 'MEME DID NOT COMPLETE WITHIN 10 MINUTES'
    else:
        fo = open( 'motif_clusters_/%04d_memeOut.txt'%(param_I_str,iii), 'r' )
        memeOut = fo.read()
        fo.close()
        #handle = StringIO.StringIO(memeOut)
        #record = motifs.parse(handle, "meme")
        #handle.close()
    record = sites = smallText = None
    try:
        record = meme.read_meme_output(memeOut, 1)[0]
        sites = pd.DataFrame(record.sites)[5].values
        smallText = str( iii ) + ' ' + str(len(clusters[iii])) + ' ' + str(record.num_sites) + \
                    ' E=' + str(record.evalue) 
    except:
        print 'COULD NOT PARSE MEME OUTPUT'
        
    return { 'sites':sites, 'record':record, 'smallText':smallText, 'memeOut':memeOut }

from cmonkeyobj import cMonkey2 as cm2
## pre-read in all cmonkey database objects for speed
cmdbs = sorted( glob.glob( '%s-out-*/cmonkey_run.db'%(organism) ) )
cms = {}
for db in cmdbs:
    print db
    try:
        cms[db] = cm2(db)
    except:
        print 'OOPS!'
        continue
#cms = { db:cm2(db) for db in cmdbs }
    
for iii in xrange( 0,len(clusters) ):
    ##df = clust_dfs.ix[iii]  # dont need for this?
    if len(clusters[iii]) < 10:
        break
    out = get_motif_cluster_sites( iii, force=False )
    if out['record'] is None:
        continue
    pdf_data = plot_motif_from_sites( out['sites'], 'pdf', out['smallText'] )
    ut.writeLines( pdf_data.split( '\n' ), 'motif_clusters_%s/%04d.pdf'%(param_I_str,iii) )
    ut.writeLines( out['memeOut'].split( '\n' ), 'motif_clusters_%s/%04d_memeOut.txt'%(param_I_str,iii) )
    print iii, 'DONE'
        
#os.popen( 'pdftk motif_clusters_%s/????.pdf cat output motif_clusters_%s/ALL.pdf'%(param_I_str,param_I_str) ).read()
try:
    os.popen( ('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=motif_clusters_%s/ALL.pdf '
               'motif_clusters_%s/????.pdf') % (param_I_str,param_I_str) ).read()
    os.popen( 'pdfnup --landscape --suffix nup --nup 5x6 motif_clusters_%s/ALL.pdf'%(param_I_str) ).read()
except:
    print 'gs or pdfnup not available!'

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

