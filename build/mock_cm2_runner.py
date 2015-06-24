#!/usr/bin/env python

import argparse

"""
A mock cm2_runner. It accepts the same input parameters and generates
the same output file, but it does nothing.
"""

DESCRIPTION = 'a mock cmonkey2 runner'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--inputfile', required=True)
    parser.add_argument('--run_num', type=int, required=True)
    parser.add_argument('--outdb', required=True)
    args = parser.parse_args()
    
    with open(args.outdb, 'w') as outfile:
        outfile.write('this is empty')
