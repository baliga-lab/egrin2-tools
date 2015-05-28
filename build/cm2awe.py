#!/usr/bin/env python

import argparse
import os
import logging
import random

import cmconfig
import ensemble
import shock
import tempfile
import json


DESCRIPTION = "cm2awe.py - prepare cluster runs for KBase AWE"

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
LOG_FILE = '/home/ubuntu/cm2awe.log'


def config2shock(targetdir):
    #outfile.write("SHOCK_URL: %s" % os.environ['SHOCK_URL'])
    #outfile.write("TOKEN: %s" % os.environ['KB_AUTH_TOKEN'])

    result = []
    for filename in os.listdir(targetdir):
        result.append(filename)
    return result

if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=LOG_FILE)
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True, help="3 letter organism code")
    parser.add_argument('--ratios', required=True, help="Path to ratios file")
    parser.add_argument('--outfile', required=True, help="file describing the Shock nodes of generated data")
    parser.add_argument('--blocks', default=None, help="Path to block definitions")
    parser.add_argument('--inclusion', default=None, help="Path to inclusion block definitions")
    parser.add_argument('--exclusion', default=None, help="Path to exclusion block definitions")
    parser.add_argument('--nruns', type=int, default=100, help="Number of cMonkey2 runs")
    args = parser.parse_args()

    targetdir = tempfile.mkdtemp(suffix='cmsplit')
    ensemble.make_ensemble_ratios(args.ratios, args.blocks, args.exclusion, args.inclusion,
                                  args.nruns, targetdir)
    cmconfig.make_config_files(1, "set1", "setfile1", args.nruns, None, targetdir)

    outinfo = config2shock(targetdir)
    with open(args.outfile, 'w') as outfile:
        outfile.write(json.dumps(outinfo))
