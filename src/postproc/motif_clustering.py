#!/usr/bin/python

import os
import bz2
import glob
import shelve
import argparse
import numpy as np
import numpy.core.defchararray as npstr
import pandas as pd
import jgraph as ig

import postproc.coding_fracs as cf

def system(cmd):
    print(cmd)
    return os.popen(cmd).read()

def main():
    parser = argparse.ArgumentParser(description="run motif clustering")
    parser.add_argument('--input_dir', required=True,
                        help="The location of tomtom results, bzip'd")
    parser.add_argument('--filter_motifs', default=True, action='store_true',
                        help="Pre-filter motifs to include; see comments")
    parser.add_argument('--plot_motifs', default=False, action='store_true',
                        help="Plot motif clusters in separate PDFs")
    parser.add_argument('--option', type=int, default=1,
                        help="Filtering option, 1 or 2; see comments")
    parser.add_argument('--mcl_I', type=float, default=2.4,
                        help='mcl I parameter value')
    parser.add_argument('--features', required=True, help='path to RSAT features file')
    args = parser.parse_args()

    # if option == 1, then we SUM all the log-pvals for multiple occurrences of the same motif pair,
    # then filter the sums to only those that are < -10 (same as done for original Halo ensemble)
    # if option == 2, then we take the LOWEST log-pval over multiple occurrences of the same motif
    # pair, with no additional filtering. Note, option (2) was used for Eco and (1) for Halo in MSB EGRIN2 paper.
    option = args.option
    print('OPTION:', option)

    # pre-filter motifs, not implemented yet. (1) remove motifs that are in coding regions (from fimo table);
    # (2) filter by motif E-value (3) filter by bicluster residual?
    pre_filter = False
    if args.filter_motifs:
        pre_filter = args.filter_motifs

    coding_fracs = total_frac = None

    if pre_filter:
        try:  # necessary only because egrin2-tools has hyphen and can't have hyphens in python module paths...
            os.symlink('egrin2-tools/postproc/coding_fracs.py', 'coding_fracs.py')
        except:
            None
        total_frac = cf.get_total_coding_rgn(args.features)

        cf_files = np.sort(np.array(glob.glob(os.path.join('*/coding_fracs.tsv.bz2'))))
        coding_fracs = []  # {}
        for f in cf_files:
            print(f)
            cff = pd.read_table(bz2.BZ2File(f), sep='\t')
            cm_run = os.path.dirname(f)  # .split('-')[2]
            cff['cm_run'] = cm_run
            if cff.shape[0] > 1:
                coding_fracs.append(cff)  # [f] = cff
        coding_fracs = pd.concat(coding_fracs, keys=None, ignore_index=True)
        # this has a hack - for some reason cluster_id in coding_fracs is %04d, trim first zero to make it %03d ...
        splitted = npstr.split(coding_fracs.motif.values.astype(str), '_')
        clust_id = np.char.mod('_%03d_', np.array([int(i[0]) for i in splitted]))  # see https://stackoverflow.com/a/28286749
        mot_id = np.char.mod('%02d', np.array([int(i[1]) for i in splitted]))
        mot_id = npstr.add(clust_id, mot_id)
        coding_fracs['motif_id'] = npstr.add(coding_fracs.cm_run.values.astype(str), mot_id)

    input_dir = 'tomtom_out'
    input_dir = args.input_dir
    files = np.sort(np.array(glob.glob(input_dir + "/*tomtom.tsv.bz2")))  # folder with the tomtom files bzip'd
    dfs = {}
    #  can pd.concat work on shelved dataframes? YES. Note protocol=2 is faster and smaller.
    dfs = shelve.open('tomtom_shelf.db', protocol=2, writeback=False)

    if len(dfs) != len(files):  # if using a shelf, once this is done once, you don't have to do it again.
        for f in files:
            gene = os.path.basename(f).split('.')[0]
            print(f, gene)
            if gene in dfs.keys():
                continue
            try:
                df = pd.read_table(bz2.BZ2File(f), sep='\t')
                print(df.shape)
                if df.shape[0] <= 0:
                    continue
                # df['gene'] = gene
                # df = pd.concat([df, pd.Series(np.repeat(gene,df.shape[0]))], 1)
                # cols = df.columns.values; cols[-1] = 'gene'; df.columns = cols
                df = df.ix[df['p-value'] <= 0.01]  # 0.1]
                print(df.shape)
                df = df.ix[df['#Query ID'] != df['Target ID']]
                print(df.shape)
                df = df.ix[df.Overlap >= 6]  # same setting as Halo run
                df = df.drop(['Query consensus', 'Target consensus'], axis=1)

                if pre_filter:
                    # add the coding fracs to the df:
                    tmp = pd.merge(df, coding_fracs, how='left', left_on='#Query ID', right_on='motif_id')
                    tmp = pd.merge(tmp, coding_fracs, how='left', left_on='Target ID', right_on='motif_id')
                    tmp = tmp.drop(['motif_x', 'cm_run_x', 'motif_id_x', 'motif_y', 'cm_run_y', 'motif_id_y'], axis=1)
                    # drop the motifs with coding fracs greater than (expected value) + (obs. stddev)/2
                    cutoff = total_frac  # + coding_fracs.coding_frac.mad() / 2
                    tmp = tmp.ix[np.logical_or(tmp.coding_frac_x.values < cutoff, tmp.coding_frac_y.values < cutoff)]
                    print(tmp.shape)
                    df = tmp; del tmp

                dfs[gene] = df

            except:
                continue

            # if f == files[20]:
            #   break

    if not os.path.isfile('motifs_tomtom.tsv.bz2'):
        if type(dfs) == dict:
            dfs2 = pd.concat(dfs, axis=0)
        else:
            dfs2 = pd.concat(dfs.values(), axis=0)
        # dfs.close()
        print(dfs2.shape)

        # incase we fail on steps below...
        dfs2.to_csv(bz2.BZ2File('motifs_tomtom.tsv.bz2', 'w'), sep='\t', index=False, header=True)

    else:
        dfs2 = pd.read_table(bz2.BZ2File('motifs_tomtom.tsv.bz2', 'r'))

    if option == 2:
        dfs2.sort('p-value', inplace=True)  # sort so lower p-values come first (these are kept by drop_duplicates)
        dfs2.drop_duplicates(['#Query ID', 'Target ID'], inplace=True)  # no, we sum up the duplicate weights below
        print(dfs2.shape)

    # all_motifs = np.sort(np.unique(np.concatenate((np.unique(dfs2['#Query ID'].values),
    #                                                   np.unique(dfs2['Target ID'].values)), 0)))
    # all_motifs_lookup = dict(itertools.izip(all_motifs.tolist(), range(len(all_motifs))))
    # inds1 = [all_motifs_lookup[i] for i in dfs2['#Query ID'].values]
    # inds2 = [all_motifs_lookup[i] for i in dfs2['Target ID'].values]

    # X = sparse.lil_matrix((len(all_motifs), len(all_motifs)))
    # X[inds1, inds2] = -np.log10(dfs['p-value'].values)

    gr = pd.DataFrame({'query': dfs2['#Query ID'].values, 'target': dfs2['Target ID'].values,
                       'weight': np.round_(-np.log10(dfs2['p-value'].values + 1e-99), 4)})
    # gr = pd.DataFrame({'query': inds1, 'target': inds2, 'weight': np.round_(-np.log10(dfs2['p-value'].values+1e-99), 4)})
    gr.to_csv('motifs_graph.tsv', sep=' ', index=False, header=False)  # igraph cannot read bzipped files (streams)
    del gr
    gr2 = ig.Graph.Read_Ncol('motifs_graph.tsv', names=True, weights=True, directed=False)
    system('pbzip2 -fv9 motifs_graph.tsv &')
    # gr2 = ig.Graph(zip(inds1, inds2), edge_attrs={'weight':np.round_(-np.log10(dfs2['p-value'].values+1e-99), 4)})
    print(gr2.ecount(), gr2.vcount())
    # gr2.write_ncol("mot_metaclustering.txt")

    # cool! see http://igraph.org/python/doc/igraph.GraphBase-class.html#simplify
    # add up the weights for duplicated edges into a single edge weight
    if option == 1:
        gr2a = gr2.simplify(multiple=True, loops=False, combine_edges=sum)
        print(gr2a.ecount(), gr2a.vcount())
        # see http://igraph.org/python/doc/tutorial/tutorial.html#selecting-vertices-and-edges
        gr2b = gr2a.es.select(weight_gt=10)  # used 10 for Halo; use less for fewer runs. returns an EdgeList
        gr2b = gr2b.subgraph()   # convert to a graph
        print(gr2b.ecount(), gr2b.vcount())
    elif option == 2:
        gr2a = gr2.simplify(multiple=True, loops=False, combine_edges=max)
        print(gr2a.ecount(), gr2a.vcount())
        gr2b = gr2a

    del gr2
    gr2b.write_ncol("mot_metaclustering.txt", weights=None)  # no weights used - same as Halo analysis which was best!
    # gr2b = ig.Graph.Read_Ncol('mot_metaclustering.txt', names=True, weights=True, directed=False)

    # now run mcl, latest version from http://www.micans.org/mcl/src/mcl-latest.tar.gz
    param_I = args.mcl_I
    cmd = './progs/mcl mot_metaclustering.txt --abc -I %.1f -v all -te 3 -S 200000' % (param_I)
    system(cmd)

    param_I_str = str(param_I).replace('.','')
    fo = open('out.mot_metaclustering.txt.I%s'%(param_I_str), 'r')
    lines = fo.readlines()
    fo.close()
    ##lines = [np.array(line.split(), np.int) for line in lines]
    lines = [np.array(line.split()) for line in lines]
    clusters = lines   ## file contains actual motif ids rather than numbers
    del lines

    # clusters = [all_motifs[i] for i in lines]
    # lines2 = ['\t'.join(i) for i in clusters]
    # fo = open('out.mot_metaclustering.txt.I45.txt', 'w')
    # fo.write('\n'.join(lines2))
    # fo.close()

    clust_lens = np.array([len(i) for i in clusters])
    print('Clusters with >= 10 motifs:', np.sum(clust_lens >= 10))
    print('Total number of motifs:', gr2a.vcount())
    print('Number of motifs in >= 10-size clusters:',
          np.sum(np.array([len(clusters[i]) for i in np.where(clust_lens >= 10)[0]])))
    print('Fraction of motifs in >= 10-size clusters:',
          float(np.sum(np.array([len(clusters[i]) for i in np.where(clust_lens >= 10)[0]]))) /
          float(gr2a.vcount()))
    del gr2a

    # Get info on alignments for each motif cluster
    dfs2.set_index('#Query ID', drop=False, inplace=True)
    clust_dfs = {}
    for i in xrange(len(clusters)):
        clust = clusters[i]
        print(i, len(clust))
        if i in clust_dfs.keys() or len(clust) < 10 or i > 500:
            continue
        df = dfs2.ix[clust]  # [np.in1d(dfs2['#Query ID'].values, clust)]
        df = df.iloc[np.in1d(df['Target ID'].values, clust)]
        df = df.sort(['p-value'])
        df = df.ix[~df.duplicated(['#Query ID', 'Target ID'])]  # remove dupes
        df = df.reset_index(drop=True)
        df['motif_clust'] = i
        print(df.shape)
        clust_dfs[i] = df

    del dfs2
    clust_dfs = pd.concat(clust_dfs, axis=0)
    print(clust_dfs.shape)
    # get coding fracs per motif cluster via:
    # clust_dfs.groupby('motif_clust').mean().coding_frac_x

    clust_dfs.to_csv(bz2.BZ2File('motif_clusts_%s.tsv.bz2'%(param_I_str), 'w'), sep='\t', index=False, header=True)
    # clust_dfs = pd.read_table(bz2.BZ2File('motif_clusts.tsv.bz2', 'r'))

    # Create the "aligned pssms" and(or?) "combined pssm" for each motif cluster
    # Not ready yet.

    if args.plot_motifs:
        # run plot_motif_clusters.py
        None
