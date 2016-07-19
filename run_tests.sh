#!/bin/bash

PYTHONPATH=. test/query_test.py && PYTHONPATH=. test/assembler_test.py && \
    PYTHONPATH=. test/util_test.py && PYTHONPATH=. test/resample_test.py
