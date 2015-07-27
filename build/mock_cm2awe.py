#!/usr/bin/env python

"""
This is a mock replacement tool for the AWE pipeline.
It takes the same input parameters and writes to the same output file
to pretend this is a actual AWE filter.

But it does nothing besides that.
"""
import argparse
import json


DESCRIPTION = 'a mock splitter tool'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="3 letter organism code")
    parser.add_argument('--ratios', required=True, help="Path to ratios file")
    parser.add_argument('--outfile', required=True, help="file describing the Shock nodes of generated data")
    parser.add_argument('--blocks', default=None, help="Path to block definitions")
    parser.add_argument('--inclusion', default=None, help="Path to inclusion block definitions")
    parser.add_argument('--exclusion', default=None, help="Path to exclusion block definitions")
    parser.add_argument('--nruns', type=int, default=100, help="Number of cMonkey2 runs")
    args = parser.parse_args()

    with open(args.outfile, 'w') as outfile:
        outfile.write(json.dumps({'this': 'is a mock test'}))
