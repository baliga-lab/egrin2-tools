#!/tools/bin/python

import optparse
import sys
import subprocess
import os
import glob
import csv
import bz2

import time

def main():

	#  Collect & check args

	op = optparse.OptionParser()
	op.add_option('-c', '--cache_dir', help='The cmpython cache directory where the genome features and info live')
	op.add_option('-f', '--features_file', help='The sequence file (with coding regions) (found in cache/<organism name>features')
	op.add_option('-o', '--organism_name', help='The organism name (e.g. eco or hal)')
	op.add_option('-i', '--input_dir', help='The cmonkey run input dir where the fimo files live (e.g. eco-out-001)')
	opt, args = op.parse_args()

	if not opt.cache_dir:
		op.error('need --cache_dir option.  Use -h for help.')
	if not opt.features_file:
		op.error('need --features_file option.  Use -h for help.')
	if not opt.organism_name:
		op.error('need --organism_name option.  Use -h for help.')
	if not opt.input_dir:
		op.error('need --input_dir option.  Use -h for help.')

	fimo_files = glob.glob(opt.input_dir + "/fimo-outs/fimo-out-*")
	print 'Number of fimo files = ', str(len(fimo_files))

	fimodict = {}

	for i in range(0,len(fimo_files)):
		fin = fimo_files[i]
		bc = os.path.basename(fin)

#		with open(fin, 'rb') as f:
		with bz2.BZ2File(fin, 'rb') as f:
			
			innerdict = {}
			cnt = 0

			for line in csv.DictReader(f, dialect='excel-tab'):
				
				strand = line['strand']
				motif = line['matched sequence']

				if strand is '+':
					stop = line['stop']
				if strand is '-':
					stop = line['start']

				innerdict[motif] = (stop, strand)

			fimodict[bc] = innerdict

		f.close()


	#  Get coding regions from features file (e.g. Escherichia_coli_K12_features)
	#  Simple CDS array of tuples containing starts and strands 
	cdsstarts = []
	cdsstops = []
	cdsstrands = []

	with open(opt.features_file, 'rb') as f:

		line = f.readline()
		while 'header' not in line:
			line = f.readline()

		for line in csv.DictReader(f, dialect='excel-tab'):
			cds = line['type']
			if cds is 'CDS':
				continue
			
			strand = line['strand']

			if strand is 'R':
				startpos = line['end_pos']
				stoppos = line['start_pos']
			else: # D or DR
				startpos = line['start_pos']
				stoppos = line['end_pos']

			cdsstarts.append(startpos)
			cdsstops.append(stoppos)
			cdsstrands.append(strand)


	f.close()

	badmotifs = {}

	for bc in fimodict: # So... there are multiple motifs in a bicluster fimo file
                print bc
		fimomotifs = fimodict[bc]
		badmotifsset = set()

		for fimomotif in fimomotifs:

			info = fimomotifs[fimomotif]
			#fimostart = int(info[0])
			fimostop = int(info[0])
			fimostrand = info[1]

			# If fimostrand is '+' look up fimo stop in cdsstarts to see if it is greater (and then check too it is less than cdsstop)
			# If fimostrand is '-' look up fimo stop in cdsstarts to see if it is less than (and also greater than cdsstop)

			for i in range(0, len(cdsstarts)):

				if (fimostrand is '+' and cdsstrands[i] is 'D') or \
						(fimostrand is '-' and cdsstrands[i] is 'R'):

					if (fimostop >= int(cdsstarts[i]) and fimostop <= int(cdsstops[i])):
						badmotifsset.add(fimomotif) # within CDS or overlapping at beginning of CDS

		badmotifs[bc] = badmotifsset
		#frac = 1.0 - (float(len(badmotifs[bc]))/float(len(fimodict[bc])))
		#print frac

	#  Print the fimo file name and the resulting fraction of goodmotifs to all motifs
        # NO THIS IS WRONG - should have separate fractions for each of the 2 motifs in each fimo file!!!

	OUT = open(os.path.join(opt.input_dir,'coding_fracs.tsv'), 'wb')

	for bc in fimodict:
		if len(fimodict[bc]) > 0:
			frac = 1.0 - (float(len(badmotifs[bc]))/float(len(fimodict[bc])))
			OUT.write(bc+'\t'+str(frac)+'\n')
		else:
			OUT.write(bc+'\tNA\n')

	OUT.close()

if __name__ == '__main__':
	start_time = float(time.time())
	main()
	print "--- %f seconds ---" % (float(time.time()) - start_time)
