![cMonkey2 Logo](https://github.com/scalefreegan/egrin2-tools/blob/master/egrin2_logo_80px.png "EGRIN2.0 Logo")

## EGRIN**2.0**-tools - Utilities and scripts for generating EGRIN**2.0** ensembles from [cMonkey2 biclustering algorithm](https://github.com/baliga-lab/cmonkey2/)

### Description

This is a collection of software tools to perform cMonkey ensemble runs and support analysis of the results.

### Documentation

A complete set of documentation for installation and running of cMonkey is on the [wiki](https://github.com/baliga-lab/cmonkey2/wiki). There are also [developer](https://groups.google.com/d/forum/cmonkey-dev) and [user](https://groups.google.com/d/forum/cmonkey-users) discussion groups. 

### System requirements

* Developed and tested with Python 2.7.6
* biopython >= 1.65
* joblib >= 0.7.1
* matplotlib >= 1.3.1
* monary >= 0.2.3
* numpy >= 1.8.2
* pandas >= 0.15.2
* scipy >= 0.13.3
* sqlite3 >= 2.6.0

for running the unit tests (optional):

* python-xmlrunner 

### Running the Unit Tests

    ./run_tests.sh
    
### Choose blocks of conditions for running a cMonkey2 ensemble

### Processing a cMonkey2 ensemble

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

### Test with 5 run *E. coli* ensemble

There is a startup script for cMonkey to run the current integrated
system

    ./cmonkey.py --organism hal --ratios example_data/hal/halo_ratios5.tsv






