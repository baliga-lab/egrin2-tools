#!/usr/bin/env python

"""cm2_runner.py - running cmonkey with input from Shock"""

import sys
import json
import os
import shock
import tempfile
import traceback
import argparse
import logging
import subprocess
import shutil


DESCRIPTION = "cm2_runner.py - run a cmonkey run using shock"

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
LOG_LEVEL = logging.DEBUG
SHOCK_URL = 'http://10.1.16.112:7078'
AUTH_TOKEN = 'un=nwportal|tokenid=b79f18d4-0ee9-11e5-b05f-123139260d4e|expiry=1465419267|client_id=nwportal|token_type=Bearer|SigningSubject=https://nexus.api.globusonline.org/goauth/keys/9269647a-0bab-11e5-bed7-123139260d4e|sig=81ad6dd81e39bde34dfda30e5689968f07211997d3fdaa02dcf585a99fc51a28e9e25bc0912fd106e8004c36c4e0ccaac22d36e3e0ad2a567f9593a3a9e77e5a2186ab0b97edae9949d6f8d2bae5741dc898ccb2785ee670666953f3176ebaef010d018d5dae24ba03c2d8ec912705d2e8349b09302741c6015c64de21f3ffe6'

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('--organism', required=True)
  parser.add_argument('--inputfile', required=True)
  parser.add_argument('--run_num', type=int, required=True)
  parser.add_argument('--outdb', required=True)
  args = parser.parse_args()

  with open(args.inputfile) as infile:
    rundata = json.load(infile)

  run_number = args.run_num
  config_node = rundata["config-%03d.ini" % run_number]
  ratio_node = rundata["ratios-%03d.tsv" % run_number]

  if 'SHOCK_URL' in os.environ:
    service_url = os.environ['SHOCK_URL']
  else:
    service_url = SHOCK_URL
  if 'KB_AUTH_TOKEN' in os.environ:
    auth_token = os.environ['KB_AUTH_TOKEN']
  else:
    auth_token = AUTH_TOKEN

  if 'LOG_DIRECTORY' in os.environ:
    logfile = os.path.join(os.environ['LOG_DIRECTORY'], 'cm2_runner-%d.log' % run_number)
    cm_logfile = os.path.join(os.environ['LOG_DIRECTORY'], 'cmonkey-%d.log' % run_number)
  else:
    logfile = None  #'/home/ubuntu/awe-logs/cm2_runner-debug.log'
    cm_logfile = '/home/ubuntu/awe-logs/cmonkey-debug.log'
  logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                      level=LOG_LEVEL, filename=logfile)


  try:
    shock_client = shock.ShockClient(service_url, auth_token)

    configfile = tempfile.NamedTemporaryFile(mode="wb", delete=False)
    ratiofile = tempfile.NamedTemporaryFile(mode='wb', delete=False)
    outdir = tempfile.mkdtemp()
    shock_client.download_node(config_node, configfile)
    shock_client.download_node(ratio_node, ratiofile)

    logging.info("Run %d", run_number)
    logging.info("shock url: %s", service_url)
    logging.info("auth token: %s", auth_token)
    logging.info("ratio node: %s, written to: %s", ratio_node, ratiofile.name)
    logging.info("config node: %s, written to: %s", config_node, configfile.name)
    logging.info("output directory: %s", outdir)
    logging.info("running cmonkey...")

    subprocess.check_output(["cmonkey.py", "--organism", args.organism,
                           "--ratios", ratiofile.name,
                           "--config", configfile.name,
                           "--out", outdir,
                           "--logfile", cm_logfile])
    logging.info("Cmonkey returned")
    # all we need now is the result databases
    # store that under a specified name
    shutil.copyfile(os.path.join(outdir, 'cmonkey_run.db'), args.outdb)
    logging.info("done.")
  except subprocess.CalledProcessError, e:
    logging.info("THE OUTPUT: %s", e.output)
    logging.exception("Error while running cmonkey runner")
    raise
  except Exception, e:
    logging.exception(e.output)
    raise

