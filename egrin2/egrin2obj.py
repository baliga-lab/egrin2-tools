import glob
import os
import numpy as np
import pandas as pd

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
        row_members = pd.concat( row_members, 0 )
        return row_members

#from egrin2.egrin2obj import EGRIN2 as eg2
#b = eg2('eco-out-')
#np.array([bb.config.getfloat('Rows','scaling_const') for bb in b.cms.values()],np.float)
