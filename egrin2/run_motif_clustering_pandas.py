#!/usr/bin/python

import bz2
import glob
import time
import math
import numpy as np
import os
from igraph import *
from pandas import *

def main():

	files = glob.glob("/Users/mharris/Documents/work/ecoli/tomtom_clust_motifs/good/*.bz2") # My test folder with the meme files bzip'd
	cnt = 1
	#pvs = []
	df = DataFrame()
	for f in files:
		print(f)
		lines = bz2.BZ2File(f).readlines()
		datamatrix = [line.split() for line in lines[1:]] # Place data in 2D array
		a = np.array(datamatrix).transpose()
		#df = DataFrame(datamatrix,columns=['Query ID','Target ID','Optimal offset','p-value','E-value','q-value','Overlap','Query consensus','Target consensus','Orientation'])
		#pv = df.iloc[:,3] # pvalues


		#dftmp = [df[:][3],df[:][7],df[:][8]]

		m1s = [a[7][i] for i in range(len(a[0])) if float(a[3][i]) <= 0.01] # Query motif, filtered by p-value
		m2s = [a[8][i] for i in range(len(a[0])) if float(a[3][i]) <= 0.01] # Target motif, filtered by p-value
		pvs = [float(el) for el in a[3] if float(el) <= 0.01 ] # Filter p-value for significant vals
		pvs = [math.log10(p) if p != 0.0 else 0.0 for p in pvs ]  # log10 unless 0, then keep as zero

		#  Add the motifs combos and pvalues to dataframe
		
		for i in range(len(m1s)):
			m1 = m1s[i]
			m2 = m2s[i]
			if m1 == m2: #  Get rid of self-self
				continue
			combo = "%s,%s" % (m1,m2)
			pv = pvs[i]

			#  The 'combo' is two motifs joined, and used as a unique column in the dataframe
			#  The p-values are added up if the combo occurs more than once - see below
			
			cols = df.columns.values
			if combo in cols:
				df[combo] = df[combo]+pv
			else:
				df[combo] = pv
			

			#df = df.append(DataFrame([pv],columns=[combo]))


		if cnt > 10:
			break
		cnt+=1

	print df.shape





if __name__ == '__main__':
	start_time = float(time.time())
	main()
	print "--- %f seconds ---" % (float(time.time()) - start_time)