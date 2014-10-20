#!/usr/bin/env python

"""This filter script takes a communities file, which is a tab separated file
with the column format

gene1 gene2 community_id community_density community_weighted_density

and returns the entries that are greater than"""

import sys

if __name__ == '__main__':
    for line in sys.stdin:
        gene1, gene2, community_id, community_density, community_weighted_density = line.split('\t')
        community_weighted_density = float(community_weighted_density)
        if community_weighted_density > 0.0:
            print line.strip()

