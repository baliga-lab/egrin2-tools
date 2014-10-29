#!/tools/bin/python
""" egrin2-fimo.py
Author:  Micheleen M. Harris
Description:  FIMO step of the EGRIN2 pipeline.  This program makes a script which can be submitted to an SGE scheduler on a cluster.

This script assumes an ensemble run under a directory structure that 
is as follows:

<dir>
  |-- <prefix>001
             |-- cmonkey_run.db
             |-- ...
  |-- <prefix>002
  |-- ...
  |-- <this is where shell scripts will be run with qsub>

Instructions:  Run resulting script (<name you choose>.sh) under <dir> as shown above.
Note:  if you are using c-shell make sure to use the --csh flag when running this program.

Example:

python egrin2-fimo.py --seqs_file cache/Escherichia_coli_K12_NC_000913.2 --user mharris --base_dir . --target_dir . --qsub_script fimo_batch_script.sh --num_cores 1  --csh --organism_name eco

Output will be in subdirectory to input directory (<prefix>001 for example will have all fimo files associated with the meme files).
"""

import optparse
import sys
import os
import glob

#  TODO:
#	-Instead of creating a fimo_script.sh for every run with simple multiple calls to fimo
#		make it such that this is more parallel (see a fimo_script.sh file under a run dir)
#	-Fix shell header to bash with an option like the '-f' in the csh shell header (maybe -profile?)
#	-Make the iteration through all MEME files to check e-value formats more efficient

#  Fimo command
FIMO_TEMPLATE = """fimo --max-stored-scores 9999999 --max-seq-length 1e8 --text --verbosity 2 %s %s | bzip2 -c > %s.bz2"""

SHELL_HEADER = """#!/bin/bash"""

# Templates for Bourne Shell
QSUB_TEMPLATE_HEADER = """#!/bin/bash

export LD_LIBRARY_PATH=/tools/lib:/tools/R-3.0.3/lib64/R/lib
export PATH=/tools/bin:${PATH}
export BATCHNUM=`printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE = """#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -M %s@systemsbiology.org
#$ -t 1-%s
#$ -cwd
#$ -pe serial %s
#$ -l mem_free=32G

%s"""


# Templates for csh

SHELL_HEADER_CSH = """#!/bin/csh -f"""

QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh -f

setenv LD_LIBRARY_PATH /tools/lib:/tools/R-3.0.3/lib64/R/lib
setenv PATH /tools/bin:${PATH}
set BATCHNUM="`printf '%03d' ${SGE_TASK_ID}`"
"""

QSUB_TEMPLATE_CSH = """#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -M %s@systemsbiology.org
#$ -t 1-%s
#$ -cwd
#$ -pe serial %s
#$ -l mem_free=32G

%s"""

def fix_meme_files(meme_files):

	newmemes = [] # The new file names - TODO:  can get rid of this by just using os to move the name to original

	#  Read memefile line by line and output line by line, replacing strange e-value formats (e.g. "E= 10.0e+003" to "E= 1.0e+004")
	for memef in meme_files:
		IN = open(memef,'rb')
		OUT = open(memef+'4fimo', 'wb')
		for line in IN:
			if 'letter-probability matrix' in line: # Note: fimo uses this matrix for input
				linespl = line.split()
				evalue = linespl[9].split('e')
				p1 = float(evalue[0]) # first part of evalue, want format to be 1.0 instead of 10.0 for e.g.
				if p1 >= 10.0: # this is something we need to change!
					p1 /= 10. # shift decimal
					p1 = round(p1,1) # round to 1 decimal place
					p2 = evalue[1] # second part of evalue, will have a + or - at beginning
					newevalue = ''
					if '+' in p2:
						p2 = p2.replace('+','')
						if p2 == '000':
							p2 = 0
							newevalue = "%se+%03d" % (str(p1),p2)
						else:
							p2 = int(p2.lstrip("0")) + 1 # now an int
							newevalue = "%se+%03d" % (str(p1),p2)
					else:
						p2 = p2.replace('-','')
						if p2 == '000':
							p2 = 0
							newevalue = "%se+%03d" % (str(p1),p2)
						else:
							p2 = int(p2.lstrip("0")) - 1 # now an int
							newevalue = "%se-%03d" % (str(p1),p2)
					linespl[9] = newevalue
					OUT.write(' '.join(linespl)+'\n')
			else:
				OUT.write(line)
		IN.close()
		OUT.close()
		os.rename((memef+'4fimo'),memef)


def main():

	#  Collect & check args

	op = optparse.OptionParser()
	op.add_option('-g', '--seqs_file', help='The sequence file (genome) (found in cache/<organism name>')
	op.add_option('-o', '--organism_name', help='KEGG name for organism (e.g. eco, hal')
	op.add_option('-u', '--user', help='User name on cluster')
	op.add_option('-i', '--base_dir', help='Cmonkey-python base directory.')
	op.add_option('-t', '--target_dir', help='The output directory name')
	op.add_option('-q', '--qsub_script', help='The script name for running fimo on cmonkey results')
	op.add_option('-n', '--num_cores', help='Number of cores to use on cluster')
	op.add_option('-s', '--csh', help='If c-shell indicate with this flag', action='store_true')
	opt, args = op.parse_args()

	if not opt.seqs_file:
		op.error('need --seqs_file option.  Use -h for help.')
	if not opt.target_dir:
		op.error('need --target_dir option.  Use -h for help.')
	if not opt.base_dir:
		op.error('need --base_dir option.  Use -h for help.')
	if not opt.num_cores:
		op.error('need --num_cores option.  Use -h for help.')
	if not opt.qsub_script:
		op.error('need --qsub_script option.  Use -h for help.')
	if not opt.organism_name:
		op.error('need --organism_name option.  Use -h for help.')
	if not opt.user:
		op.error('need --user option.  Use -h for help.')
	if not opt.num_cores:
		op.error('need --num_cores option.  Use -h for help.')

	if opt.csh:
		header = QSUB_TEMPLATE_HEADER_CSH
		template = QSUB_TEMPLATE_CSH
		shellheader = SHELL_HEADER_CSH
	else:
		header = QSUB_TEMPLATE_HEADER
		template = QSUB_TEMPLATE
		shellheader = SHELL_HEADER


	#  Fix genome seqs file for fimo (to upper and add header), ouput to new file, put in input dir (might be better way...)
	seqsfile_in = open(opt.seqs_file, 'rb')
	lines = seqsfile_in.readlines()
	genomeseq = lines[0].upper() # to upper may not be necessary

	seqsfile_out = os.path.join(opt.seqs_file)+'4fimo'
	out = open(seqsfile_out, 'wb')
	out.write(">"+os.path.basename(opt.seqs_file)+"\n") # The short little header needed by fimo
	out.write(genomeseq)
	out.close()

	#  Create a dict of run output dirs with array of meme file names
	org_out_dirs = glob.glob(os.path.join(opt.base_dir,"%s-out-*" % opt.organism_name))

	out_dir_dict = {}

	for org_dir in org_out_dirs:
		#  Find MEME files
		meme_files = glob.glob(os.path.join(org_dir, "meme-out-*"))
		fix_meme_files(meme_files) # Fixes e-value formats (iterates thru every line of MEME file - TODO: make more efficient)
		out_dir_dict[org_dir] = meme_files

	#  We will make a subscripts for each run directory to call fimo on all meme files within that directory
	sub_scripts = {}
	for org_dir in org_out_dirs:
		meme_files = out_dir_dict[org_dir]
		fimo_cmd = ""

		#  We will make a subscript for each run directory to call fimo on all meme files within that directory
		for meme in meme_files:
			num = os.path.basename(meme).split('-')[3] # Get the run number
			num = num[1:] # remove leading 0
			fimo_cmd += '\n' + FIMO_TEMPLATE % (meme, seqsfile_out, os.path.join(org_dir, "fimo-out-%s" % num)) + '\n'

		sub_scripts[org_dir] = fimo_cmd

	for org_dir in org_out_dirs:
		outfile_name = os.path.join(org_dir, 'fimo_script.sh')
		with open(outfile_name, 'w') as outfile:
			subscript_template = shellheader + sub_scripts[org_dir]
			outfile.write(shellheader)
			outfile.write(sub_scripts[org_dir])
			#rendered = Template(subscript_template).render()
			#outfile.write(rendered)
			#st = os.stat(outfile)
			#os.chmod(subscript, st.st_mode | stat.S_IEXEC)
		outfile.close()
		#st = os.stat(outfile)
		#os.chmod(subscript, st.st_mode | stat.S_IEXEC)
		os.chmod(outfile_name,0744)



	#  Create master script
	with open(os.path.join(opt.base_dir, opt.qsub_script), 'w') as outfile:
		if opt.user is not None:
			login = opt.user
		else:
			os.getlogin()

		num_runs = len(org_out_dirs)

		subscript = "%s-out-${BATCHNUM}/fimo_script.sh" % (opt.organism_name)


		outfile.write(header)
		outfile.write(template % (login, 
			login,
			num_runs,
			opt.num_cores, 
			subscript))

	outfile.close()


if __name__ == '__main__':
	main()
	














