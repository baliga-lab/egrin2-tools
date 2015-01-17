![cMonkey2 Logo](https://github.com/baliga-lab/cmonkey2/blob/master/graphics/cmonkey2_logo_80px.png "cMonkey2 Logo")

## EGRIN**2**-tools - Utilities and scripts for generating EGRIN<sub>2</sub> ensembles using [cMonkey2-python biclustering algorithm](https://github.com/baliga-lab/cmonkey2/)

### Description

This is a collection of software tools to perform cMonkey ensemble runs and support analysis of the results.

### Documentation

A complete set of documentation for installation and running of cMonkey is on the [wiki](https://github.com/baliga-lab/cmonkey2/wiki). There are also [developer](https://groups.google.com/d/forum/cmonkey-dev) and [user](https://groups.google.com/d/forum/cmonkey-users) discussion groups. 

### System requirements

* Developed and tested with Python 2.7.2
* scipy >= 0.9.0
* numpy >= 1.6.0
* biopython >= 1.63
* MySQLdb >= 1.2.3
* BeautifulSoup >= 3.2.0
* R >= 2.14.1
* rpy2 >= 2.2.1
* MEME 4.3.0 or >= 4.8.1
* csh (for running MEME)
for the human setup, Weeder 1.4.2 is needed

for running the unit tests (optional):

* python-xmlrunner 

for running the monitoring application (optional):

* CherryPy 3
* Jinja2
* python-routes

### Running the Unit Tests

    ./run_tests.sh

### Running cmonkey2

** Usage

The tools are intended to be run in the following order

1. Ensemble generation (egrin2/ensemble.py)
2. Run ensemble on cluster
3. Tomtom job generation
4. Run tomtom on cluster
5. Run corems

In general, you should be able to run cmonkey2 on microbial gene
expressions with

    ./cmonkey.py --organism <organism-code> --ratios <tab separated file of gene expressions>

The file can be either in your file system or a web URL.

After the program was started, a log file will be written in cmonkey.log. You
can see all available options with

    ./cmonkey.py --help

### Test Run with Halobacterium Salinarum

There is a startup script for cMonkey to run the current integrated
system

    ./cmonkey.py --organism hal --ratios example_data/hal/halo_ratios5.tsv

### Start the python based monitoring application

    python cmviewer/main.py





