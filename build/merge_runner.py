#!/usr/bin/env python
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--ncbi_code', required=True)
    parser.add_argument('--ratios', required=True)
    args = parser.parse_args()
    
