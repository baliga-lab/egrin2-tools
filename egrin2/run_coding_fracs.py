import optparse
import os
import glob
import bz2

import pandas as pd
import numpy as np

def main():
    #  Collect & check args
    op = optparse.OptionParser()
    op.add_option('-c', '--cache_dir', help='The cmpython cache directory where the genome features and info live')
    op.add_option('-f', '--features_file', help='The features file (with coding regions) (found in cache/<organism name>features')
    op.add_option('-o', '--organism_name', help='The organism name (e.g. eco or hal)')
    op.add_option('-i', '--input_dir', help='The cmonkey run input dir where the fimo files live (e.g. eco-out-001)')
    opt, args = op.parse_args()

    if not opt.cache_dir:
        op.error('need --cache_dir option.  Use -h for help.')
    if not opt.features_file:
        op.error('need --features_file option.  Use -h for help.')
    if not opt.organism_name:
        op.error('need --organism_name option.  Use -h for help.')
    if not opt.input_dir:
        op.error('need --input_dir option.  Use -h for help.')

    fimo_files = np.sort( np.array( glob.glob(os.path.join(opt.input_dir, "fimo-outs/fimo-out-*")) ) )
    print 'Number of fimo files = ', str(len(fimo_files))

    #  Get coding regions from features file (e.g. Escherichia_coli_K12_features)
    f = open(opt.features_file, 'r')
    skip = 1
    line = f.readline()
    while 'header' not in line:
        line = f.readline()
        skip += 1
    f.close()

    features = pd.read_table(opt.features_file,engine='python',skiprows=skip)
    features = features[ features.type != 'SEQ_END' ]

    for f in fimo_files:
        print f
        fimo = pd.read_table( bz2.BZ2File(f) )
        is_bad = np.zeros( fimo.shape[0], dtype=bool )
        for i in xrange(fimo.shape[0]):
            row = fimo.ix[i]
            hits = features[ (features.start_pos <= row.start) & (features.end_pos >= row.start) ]
            if hits.shape[0] <= 0:
                hits = features[ (features.start_pos <= row.stop) & (features.end_pos >= row.stop) ]
            if hits.shape[0] > 0:
                is_bad[i] = True

        fimo['in_coding_rgn'] = is_bad
        fimo.to_csv( bz2.BZ2File( f, 'w' ), sep='\t', index=False )

if __name__ == '__main__':
    main()
