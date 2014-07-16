#!/usr/bin/env python

EVALUE_CUTOFF = 100
RESID_CUTOFF  = 0.8
DIST_METHOD   = "ed"
Q_THRESHOLD   = 0.5
MIN_OVERLAP   = 4
Q_PSEUDO      = 0
T_PSEUDO      = 0

CMD_TEMPLATE = "%s/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"

if __name__ == '__main__':
    pass
