import glob
import os
import numpy as np

from cmonkeyobj import cMonkey2 as cm2

class EGRIN2:
    cms = {} ## f: cm2(f) for f in files }
    prefix = 'eco-out-'
    files = None

    def __init__( self, prefix ):
        self.prefix = prefix
        self.files = np.sort( np.array( glob.glob( prefix + "???/cmonkey_run.db" ) ) ) # get all cmonkey_run.db files

        for f in self.files:
            print f
            try:
                b = cm2(f)
                self.cms[int(os.path.dirname(f).replace(self.prefix,''))] = b
            except:
                print 'ERROR'



#from egrin2.egrin2obj import EGRIN2 as eg2
#b = eg2('eco-out-')
