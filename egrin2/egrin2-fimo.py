#!/tools/bin/python

# egrin2-fimo.py
# Author:  Micheleen M. Harris
# Description:  Reimplemtation of egrin2_fimo.R in python

import optparse
import sys
import subprocess
import os
import glob

#  TODO:  Logs, JSON output, bzip fimo output folder (if need be)

def main():

	#  Collect & check args

	op = optparse.OptionParser()
	op.add_option('-c', '--cache_dir', help='The cmpython cache directory where the genome features and info live.')
	op.add_option('-g', '--seqs_file', help='The sequence file (genome) (found in cache/<organism name>')
	op.add_option('-i', '--input_dir', help='Cmonkey-python output directory.  Input for FIMO.  Must have MEME results.')
	op.add_option('-o', '--output_dir', help='The output directory name')
	op.add_option('-p', '--progs_dir', help='The location of fimo (MEME suite bin dir)')
	op.add_option('-w', '--overwrite', help='true or false if you would like to overwrite previous fimo results')
	opt, args = op.parse_args()

	if not opt.cache_dir:
		op.error('need --cache_dir option')
	if not opt.seqs_file:
		op.error('need --seqs_file option')
	if not opt.output_dir:
		op.error('need --output_dir option')
	if not opt.input_dir:
		op.error('need --input_dir option')
	if not opt.progs_dir:
		op.error('need --progs_dir option')
	if not opt.overwrite:
		op.error('need --overwrite option')

	#  Get all meme output

	meme_files = glob.glob(os.path.join(opt.input_dir, "meme-out-*"))
	print 'Number of meme files (biclusters) = ', str(len(meme_files))

	if not os.path.exists(opt.output_dir):
		os.mkdir(opt.output_dir)

	#  Fix genome seqs file for fimo (to upper and add header), ouput to new file, put in cache dir (might be better way...)
	seqsfile_in = open(os.path.join(opt.cache_dir, opt.seqs_file), 'rb')
	lines = seqsfile_in.readlines()
	genomeseq = lines[0].upper() # to upper may not be necessary
	seqsfile_out = os.path.basename(opt.seqs_file)+'4fimo'
	out = open(os.path.join(opt.cache_dir, seqsfile_out), 'wb')
	out.write(">"+os.path.basename(opt.seqs_file)+"\n") # The short little header needed by fimo
	out.write(genomeseq)
	out.close()

	if 'true' in opt.overwrite:

		#  A dict for bicluster fimo results - for JSON (not doing presently)
		#resultsdict = {}

		#  Run FIMO on all meme files - Should I check for latest MEME suite???

		for f in meme_files:
			clust = os.path.basename(f).split('-')[3]
			outfile = os.path.join(opt.output_dir, 'fimo-out-'+clust)
			cmd = os.path.join(opt.progs_dir, "fimo") + " --max-stored-scores 9999999 --max-seq-length 1e8 --text --verbosity 2 %s %s > %s"  % (f, os.path.join(opt.cache_dir, seqsfile_out), outfile)
			subprocess.call(cmd, shell=True) #  NB:  At this point writing out fimo results to a file 'fimo-out-'+cluster number

			#  Add this fimo out to a dict obj by bicluster (name according to run dir) and jsonify after this loop (not doing presently)

			#tmpdict = readFimoOut(outfile)
		"""
		#  Remove any old bzip'd folder
		cmd = "rm -f "+os.path.join(opt.output_dir, "fimo-out.tar.bz2")
		subprocess.call(cmd,shell=True)

		#  Tar and Bzip2 fimo folder and store in 'fimo_out' folder

		cmd = "tar -cjf " + os.path.join('fimo_out', "fimo-out.tar.bz2") + " fimo_out"
		subprocess.call(cmd, shell=True)
		cmd = "mv "+os.path.join(opt.input_dir, "fimo-out.tar.bz2") + " " + opt.output_dir
		subprocess.call(cmd, shell=True)
		"""

if __name__ == '__main__':
	main()

	













