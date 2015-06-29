#!/usr/bin/env python
import argparse
import logging
import os

"""
merge_runner.py

merges all the cmonkey runs in previous steps and builds a database of Corems
(co-regulated modules).

It will then generate input for a set of resampling runs
"""
DESCRIPTION = "merge_runner.py - merge multiple cmonkey results"

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--ratios', required=True)
    parser.add_argument('result_dbs', nargs='*')
    args = parser.parse_args()

    if 'LOG_DIRECTORY' in os.environ:
        logfile = os.path.join(os.environ['LOG_DIRECTORY'], 'merge_runner.log')
    else:
        logfile = None
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=LOG_LEVEL, filename=logfile)

    logging.info("Running merge...")
    logging.info(args.result_dbs)

