import glob
import os
import numpy as np
import pandas as pd
import bz2
import sqlite3

from cmonkeyobj import cMonkey2 as cm2

class EGRIN2:
    cms = {} ## f: cm2(f) for f in files }
    prefix = 'eco-out-'
    files = None
    tables = None

    def __init__( self, prefix, nruns=None ):
        self.prefix = prefix
        self.files = np.sort( np.array( glob.glob( prefix + "???/cmonkey_run.db" ) ) ) # get all cmonkey_run.db files

        for f in self.files:
            print f
            try:
                b = cm2(f)
                self.cms[int(os.path.dirname(f).replace(self.prefix,''))] = b
            except:
                print 'ERROR: %s' % (f)

            if nruns is not None and len(self.cms) >= nruns:
                break

        self.tables = { tname: self.concat_tables(tname) for tname in self.cms[1].tables.keys() }

    def concat_tables( self, tname ):
        row_members = { k: self.cms[k].tables[tname] for k in self.cms.keys() }
        for i in row_members.keys():
            tmp = row_members[i]
            tmp = pd.concat( [ tmp, pd.Series(np.repeat(i,tmp.shape[0])) ], 1 )
            cols = tmp.columns.values; cols[-1] = 'RUNID'; tmp.columns = cols
            row_members[i] = tmp
        row_members = pd.concat( row_members, 0, ignore_index=True, copy=False )
        return row_members

    def make_dbfiles( self ):
        ## read in fimo tables and store as one big sqlite database, one for each run
        for f in self.files:
            print(f)
            fimo_outs = np.sort( np.array( glob.glob(os.path.dirname(f)+'/fimo-outs/fimo-out-*bz2') ) )
            if len(fimo_outs) <= 0:
                continue
            dfs = {}
            for ff in fimo_outs:
                print '    ' + ff
                try:
                    df = pd.read_table( bz2.BZ2File(ff), sep='\t' )
                    clust_num = int(os.path.basename(ff).split('-')[2].replace('.bz2',''))
                    df = pd.concat( [ df, pd.Series(np.repeat(clust_num,df.shape[0])), \
                                      pd.Series(np.repeat(os.path.dirname(f),df.shape[0])) ], 1 )
                    cols = df.columns.values; cols[-2] = 'cluster'; cols[-1] = 'run_id'; df.columns = cols
                    dfs[ clust_num ] = df ##.append( df )
                except:
                    print 'ERROR'
            dfs = pd.concat( dfs, 0, ignore_index=True, copy=False )
            conn = sqlite3.connect( os.path.dirname(f)+'/fimo-outs.db' )
            dfs.to_sql( 'fimo_out', conn, index=False, if_exists='replace' )
            conn.execute( 'create index idx on fimo_out(start,stop,"p-value")' )
            tmp = pd.read_sql( 'select * from fimo_out limit 10', conn )
            print tmp
            conn.close()

    def get_motifs_at_posn( self, posn, pvalue=1e-5, seqid=None ):
        ## Use the fimo databases created with make_dbfiles() to select the fimo hits within the given 
        ## coords on the given genome sequence.
        query = 'select * from fimo_out where start <= ? and stop >= ? and "p-value" <= ?'
        params = [posn, posn, pvalue]
        if seqid is not None:
            query = 'select * from fimo_out where start <= ? and stop >= ? and "p-value" <= ? and "sequence name" == seqid'
            params = [posn, posn, pvalue, seqid]

        hits = {}
        for f in self.files:
            print(f)
            try:
                conn = sqlite3.connect( os.path.dirname(f)+'/fimo-outs.db' )
                tmp = pd.read_sql( query, conn, params=params )
                conn.close()
                hits[os.path.dirname(f)] = tmp
            except:
                None
        hits = pd.concat( hits, 0, ignore_index=True, copy=False )

#from egrin2.egrin2obj import EGRIN2 as eg2
#b = eg2('eco-out-')
#np.array([bb.config.getfloat('Rows','scaling_const') for bb in b.cms.values()],np.float)
## Here's how to get nrows for all clusters in the ensemble and plot a histogram:
#tmp=b.tables['row_members'].groupby(['cluster','RUNID']).size()
#pd.DataFrame({'size':tmp}).plot(kind='hist', bins=20)
## And motifs within a location and pvalue. First (only once, do the first line):
#b.make_dbfiles()
#b.get_motifs_at_posn( 10000, 1e-5 )
