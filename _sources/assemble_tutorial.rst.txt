Assemble Tutorial
=================

**Important!!! This tutorial assumes you have access to a complete cMonkey2 ensemble.**

In a nutshell
-------------

The ASSEMBLE scripts transfer and compile individual cMonkey2 SQLite databases into an integrated MongoDB database.

In addition, they perform several post-processing steps, including: detection of gene regulatory elements (GREs) by comparing individual bicluster motifs with TOMTOM and clustering with MCL, genome-wide scanning of motifs with FIMO, and detection of co-regulated modules or **corems** using link-community detection.

Requirements
------------

  * MongoDB >= 2.4.9
  * compiled C++ scripts for corem detection, available here

IMPORTANT: This tutorial currently assumes that ``TOMTOM``, ``MCL`` and ``FIMO`` have already been run.

A single GRE definition file is read from, eg:

.. highlight:: none

::

   /ensemble-head-dir
       /out.mot_metaclustering.txt.I45.txt

FIMO scans are read from each run sub-directory, eg:

.. highlight:: none

::

   /ensemble-head-dir
       /org-out-xxx
           /fimo-outs
               /fimo-out-xxxx.bz2

Optional:

  * ``row_annot``: tab-delimited row (gene) annotations. Will be downloaded from MicrobesOnline automatically using ``--ncbi_code`` if undefined
  * ``col_annot``: tab-delimited column (condition) annotations.

The format for these files will be described in detail below.

Additionally, the Python modules described on the Home page are required to run these scripts.

Scripts
-------

  * ``assembler.py``: The control function for ASSEMBLE scripts.
  * ``makeCorems.py``: Identifies corems using C++ scripts compiled above
  * ``resample_QSub.py``: Generates QSub script for submission of resamples to cluster
  * ``sql2mongoDB.py``: Merges individual cMonkey SQLite dbs and post-processing data into MongoDB


.. highlight:: none

::

   %run ..//assembler.py -h

   usage: assembler.py [-h] --organism ORGANISM --ratios RATIOS --targetdir TARGETDIR --ncbi_code NCBI_CODE [--cores CORES] [--ensembledir ENSEMBLEDIR] [--col_annot COL_ANNOT] [--host HOST] [--port PORT] [--prefix PREFIX] [--row_annot ROW_ANNOT] [--row_annot_matchCol ROW_ANNOT_MATCHCOL] [--gre2motif GRE2MOTIF] [--db DB] [--genome_annot GENOME_ANNOT] [--backbone_pval BACKBONE_PVAL] [--link_comm_score LINK_COMM_SCORE] [--link_comm_increment LINK_COMM_INCREMENT] [--link_comm_density_score LINK_COMM_DENSITY_SCORE] [--corem_size_threshold COREM_SIZE_THRESHOLD] [--n_resamples N_RESAMPLES] [--cluster CLUSTER] [--finish_only FINISH_ONLY] [--user USER]

   assemble.py - prepare cluster runs

   optional arguments:
     -h, --help            show this help message and exit
     --organism ORGANISM   3 letter organism code
     --ratios RATIOS       Path to ratios file. Should be 'raw' (normalized) ratios, not the standardized ratios used by cMonkey
     --targetdir TARGETDIR Storage path for MongoDB and corem data
     --ncbi_code NCBI_CODE NCBI organism code
     --cores CORES Number  local cores to use for corem C++ scripts
     --ensembledir ENSEMBLEDIR
                           Path to ensemble runs. Default: cwd
     --col_annot COL_ANNOT Tab-delimited file with experiment annotations
     --host HOST           MongoDB host. Default 'localhost'
     --port PORT           MongoDB port
     --prefix PREFIX       Ensemble run prefix. Default: organism-out-
     --row_annot ROW_ANNOT Optional row (gene) annotation tab-delimited file. If not specified, annotations will be downloaded from MicrobesOnline using --ncbi_code.
     --row_annot_matchCol ROW_ANNOT_MATCHCOL
                           Name of column in row_annot that matches row names in ratios file.
     --gre2motif GRE2MOTIF Motif->GRE clustering file
     --db DB               Optional ensemble MongoDB database name
     --genome_annot GENOME_ANNOT
                           Optional genome annotation file. Automatically downloaded from MicrobesOnline using --ncbi_code
     --backbone_pval BACKBONE_PVAL
                           Significance pvalue for gene-gene backbone. Default = 0.05.
     --link_comm_score LINK_COMM_SCORE
                           Scoring metric for link communities
     --link_comm_increment LINK_COMM_INCREMENT
                           Height increment for cutting agglomerative clustering of link communities
     --link_comm_density_score LINK_COMM_DENSITY_SCORE
                           Density score for evaluating link communities
     --corem_size_threshold COREM_SIZE_THRESHOLD Defines minimum corem size. Default = 3.
     --n_resamples N_RESAMPLES
                           Number resamples to compute for corem condition assignment. Default = 10,000
     --cluster CLUSTER     Run re-samples on cluster? Boolean.
     --finish_only FINISH_ONLY
                           Finish corems only. In case session gets dropped
     --user USER Cluster   user name


ASSEMBLE an EGRIN 2.0 ensemble
------------------------------

In this tutorial we will see how you would ASSEMBLE an *Escherichia coli* EGRIN 2.0 ensemble using several example files and a couple of cMonkey2 runs, which we provide here.

STEP 1: Generate optional input files
-------------------------------------

First, let's explore the optional annotation files. Providing annotations for genes and conditions is a great way to enrich your analysis of the ensemble. You can get a better idea for the utility of these metainformation by following the advanced mining tutorial

``row_annot``
~~~~~~~~~~~~~

As noted above, the ``row_annot`` file will be downloaded automatically from MicrobesOnline if a custom annotation is not provided. If you provide your own row_annot file, however, you will also need to specificy ``--row_annot_matchCol``, which is the name of the column in your annotation file that matches the gene name used by cMonkey2 (i.e. the row names in your ratios file).

The row annotation file should look like the annotation file supplied by MicrobesOnline, where each row specifies a gene and each of the columns specifies some information about that gene. Again, you must ensure that at least one of the columns contains gene names that match the gene names in the ratios file used by cMonkey2, in the case of MicrobesOnline, it is the ``sysName`` column.

Here is an example annotation file for *E. coli* direct from MicrobesOnline, the file itself is available here.

.. figure:: _static/assemble/row_annot.png
            :alt: Example row_annot file

``col_annot``
~~~~~~~~~~~~~

The col_annot file provides metainformation about each experiment. Like the row_annot file, these annotations are optional, but they can be valuable for making sense of ensemble predicitions.

Please note that the file format is different here. Each row contains a particular experimental meta-annotation followed by several required descriptions: (1) experiment_name, (2) feature_name, (3) value, (4) feature_units, (5) feature_type. The experiment_name column should match the experiment name in the ratios file.

The `col_annot` file should look like the tab-delimited file depicted below. You can download an *E. coli* `col_annot` file to use as a template here

.. figure:: _static/assemble/col_annot.png
            :alt: Example col_annot file


STEP 2: Run ``assembler.py``
----------------------------

**Important!!! TOMTOM, MCL, and FIMO should be run prior to assembly. Otherwise the ensemble will not contain GREs or motif scans.**

To run the assembler, you must supply several files as well as specify where you would like to host the MongoDB database. If you do not supply a host for the MongoDB databse, it will be stored on ``localhost``. At a minimum, you should supply:

  * ``--organism``: 3 letter organism code
  * ``--ratios``: Path to ratios file. Should be 'raw' (normalized) ratios, not the standardized ratios used by cMonkey
  * ``--targetdir``: Storage path for MongoDB and corem data
  * ``--ncbi_code``: NCBI organism code
  * ``--ensembledir``: Path to ensemble runs. Default: cwd

Assuming your terminal current working directory is within the E. coli ensemble, with the cMonkey2 runs located in the subfolder ``./eco-ens-m3d/``, you can call the assembler as follows:

.. highlight:: none

::

   python assembler.py --organism eco --ratios ./ratios_eco_m3d.tsv.gz --targetdir ./ --ncbi_code 511145 --ensembledir ./eco-ens-m3d/ --col_annot ./E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0

We have included a small test ensemble (need to find a place to host this!!!). In this call we set ``--n_resamples`` to 0. This means that you will not have to perform condition resampling for corems, which currently requires access to an SGE cluster.

The MongoDB Database Schema is available here

STEP 3: Run ``col_resample`` on SGE cluster
-------------------------------------------

After loading all of the data from cMonkey2 into MongoDB and detecting corems, the ``assembler.py`` script will pause so that a computationally intensive procedure can be run on an SGE cluster.

This step assigns conditions to the corems. In short, it computes brute force resamples of the gene expression data (default: 10,000 resamples) to determine the conditions in which genes from each corem are tightly co-expressed. For more information about this step, please refer to the `Supplementary Information of our 2014 paper <http://msb.embopress.org/highwire/filestream/49752/field_highwire_a_enclosures/0/supplementary-material.inline-supplementary-material-10.pdf?download=true>`_.

The script generates an output directory containing QSub scripts and a master control script called ``resample.sh``. You should transfer this folder as well as the ``resample.py`` script located in the assemble folder of the egrin2-tools repository to the cluster to compute the resamples.

When ``assembler.py`` pauses for this step, you will see the following message displayed on the screen.

.. highlight:: none

::

   Output Qsub scripts to targetdir

   Transfer these documents to the cluster. Run 'resample.sh' with resample.py in your working directory to compute all resamples.

   Once this is done, return here to finish processing corems.

   Please type: 'Done' to continue


Be sure that the SGE cluster on which the resampling is run has access to the MongoDB host, as the script directly writes the output of the resampling to this database. Once the resampling is finsihed (all resample database entries present), you can return to this prompt and type Done to finish the assembly

**Alternatively:** if for some reason your session terminates prematurely or for any other reason there is an interruption at this step, you can start ``assembler.py`` with the flag ``--finish_only True``. You must also include all of the original parameters to the ``assembler.py`` script.

Importantly, this still assumes that all of the information from the resampling procedure is contained in the MongoDB database.

For example, returning to our original ``assembler.py`` call, we would type:

.. highlight:: none

::

   python assembler.py --organism eco --ratios ./ratios_eco_m3d.tsv.gz --targetdir ./ --ncbi_code 511145 --ensembledir ./eco-ens-m3d/ --col_annot ./E_coli_v4_Build_6.experiment_feature_descriptions.tsv.gz --n_resamples 0 --finish_only True

To finish the assembly.

STEP 4: Finish assembly
-----------------------

Having finished assembly, the ``assembler.py`` script will dump the MongoDB database to BSON, which are subsequently compressed in the supplied ``--targetdir``. This file can be uncompressed and read directly into MongoDB using the ``mongorestore`` command, e.g.:

.. highlight:: none

::

   mongorestore eco_db
