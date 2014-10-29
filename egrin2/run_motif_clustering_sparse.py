#!/usr/bin/python

import bz2
import glob
import time
import math
import numpy as np
import os
from scipy.sparse import *
from scipy import *
#import shelve
from igraph import *

def main():


	files = glob.glob("/Users/mharris/Documents/work/ecoli/tomtom_clust_motifs/good/*.bz2") # My test folder with the meme files bzip'd
	orgname = 'eco'

	motifs1 = []
	motifs2 = []
	pvs = []
	rowidxs = []
	colidxs = []

	idxrow = 0
	idxcol = 0



	cnt = 0

	for f in files:

		print(f)

		lines = bz2.BZ2File(f).readlines()
		datamatrix = [line.split() for line in lines[1:]] # Place data in 2D array
		a = np.array(datamatrix).transpose() # Get it so that for instance row 1 is Query ID, row 2 is Target ID, etc.
		m1 = [a[7][i] for i in range(len(a[0])) if float(a[3][i]) <= 0.01] # Query motif, filtered by p-value
		m2 = [a[8][i] for i in range(len(a[0])) if float(a[3][i]) <= 0.01] # Target motif, filtered by p-value
		motifs1.extend(m1)
		motifs2.extend(m2)
		pv = [float(el) for el in a[3] if float(el) <= 0.01 ] # Filter p-value for significant vals
		pv = [math.log10(p) if p != 0.0 else 0.0 for p in pv ]  # log10 unless 0, then keep as zero
		pvs.extend(pv)

		if cnt > 10:
			break

		cnt+=1

	print "Creating mapping..."

	#  Create a mapping of motif pairs to the pvalues - Room for improvement, maybe a db instead
	motifs2pvals = {}
	for i in range(len(pvs)):
		motifs2pvals[(motifs1[i],motifs2[i])] = pvs[i]

	#  Combine motif lists
	motifs = motifs1
	motifs.extend(motifs2)
	motifs = list(set(motifs))


	x = csc_matrix( (len(motifs),len(motifs)), dtype='d' ) # empty matrix
	x = x.tolil() # Convert to lil_matrix for quicker access to matrix elements

	print "Assigning p-values to sparse matrix..."
	start_time = float(time.time())

	for i in range(len(motifs)): # Room for improvement
		m1 = motifs[i]
		for j in range(len(motifs)):
			m2 = motifs[j]
			if m1 != m2: # Filter out self-self
				if (m1,m2) in motifs2pvals:
					pv = motifs2pvals[(m1,m2)]
					#print pv
					x[i,j] = x[i,j] + pv

	print x.get_shape()

	print "--- %f seconds ---" % (float(time.time()) - start_time)

	#  Igraph
	g = Graph()

	#  Vertices
	g.add_vertices(len(motifs))

	print "Adding edges to igraph object..."

	start_time = float(time.time())

	
	# Edges
	edges = []
	#  If there is a connection with log pvalue < -10.0 add edge
	for i in range(len(motifs)):
		edges.extend([(i,j) for j in range(len(motifs)) if x[i,j] < -10.0])

	# Simplify the graph by removing self-loops and/or multiple edges
	print "Simplifing graph...."
	g = g.simplify()

	print "--- %f seconds ---" % (float(time.time()) - start_time)






	"""

	OUT = open('mcl_input.txt', 'wb')

	print "Writing to file...."

	for i in range(len(motifs)):
		row = x.getrow(i).toarray()[0]
		row = [row[i] for i in range(len(row)) if row[i] < -10.0] # Filter by cumulative p-values (pvalues were added together for multipe motif pairs occuring in different files)
		for j in range(len(row)):
			#if row[j] < -10.0:
				#print motifs[i]
				#print motifs[j]
			OUT.write('%s\t%s\n'%(motifs[i], motifs[j]))

	OUT.close()
	"""





if __name__ == '__main__':
	start_time = float(time.time())
	main()
	print "--- %f seconds ---" % (float(time.time()) - start_time)