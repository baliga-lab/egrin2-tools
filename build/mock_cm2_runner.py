#!/usr/bin/env python

import argparse
import os
import traceback
import shutil

"""
A mock cm2_runner. It accepts the same input parameters and generates
the same output file, but it does nothing.
"""

DESCRIPTION = 'a mock cmonkey2 runner'

LOGPATH = '/home/ubuntu/awe-logs'

if __name__ == '__main__':
    logfile = open(os.path.join(LOGPATH, 'mock_runner'), 'a')
    try:
      logfile.write('starting runner\n')
      parser = argparse.ArgumentParser(description=DESCRIPTION)
      parser.add_argument('--organism', required=True)
      parser.add_argument('--inputfile', required=True)
      parser.add_argument('--run_num', type=int, required=True)
      parser.add_argument('--outdb', required=True)
      args = parser.parse_args()

      logfile.write("copying file '%s'\n" % args.outdb)
      shutil.copyfile(os.path.join('/home/ubuntu', args.outdb), args.outdb)
    except:
        logfile.write('exception\n')
        traceback.print_exc(logfile)
