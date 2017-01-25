#!/bin/bash

PYTHONPATH=`pwd` test/query_test.py $@ && PYTHONPATH=`pwd` test/assembler_test.py $@ && \
    PYTHONPATH=`pwd` test/util_test.py $@ && PYTHONPATH=`pwd` test/resample_test.py $@
