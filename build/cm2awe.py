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


def config2shock(targetdir):
    service_url = os.environ['SHOCK_URL']
    auth_token = os.environ['KB_AUTH_TOKEN']
    logging.info("SHOCK SERVICE: %s", service_url)
    logging.info("AUTH TOKEN: %s", auth_token)

    shock_client = shock.ShockClient(service_url, auth_token)
    result = {}
    for filename in os.listdir(targetdir):
        #result.append(filename)
        path = os.path.join(targetdir, filename)
        try:
            logging.info("uploading input file '%s' ...", path)
            shock_result = shock_client.upload_file(path)
            shock_node_id = shock_result['data']['id']
            result[filename] = shock_node_id
        except:
            logging.error("error uploading file '%s'", path)
            traceback.print_exc()

    return result

if __name__ == '__main__':
    if 'LOG_DIRECTORY' in os.environ:
        logfile = os.path.join(os.environ['LOG_DIRECTORY'], 'cm2awe.log')
    else:
        logfile = None
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=logfile)
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
    try:
        ensemble.make_ensemble_ratios(args.ratios, args.blocks, args.exclusion, args.inclusion,
                                    args.nruns, targetdir)
        cmconfig.make_config_files(1, "set1", "setfile1", args.nruns, None, targetdir)

        outinfo = config2shock(targetdir)
        with open(args.outfile, 'w') as outfile:
            outfile.write(json.dumps(outinfo))
    except:
        logging.exception("Exception while running")
        raise
