Postprocessing
==============

Introduction
------------

After the individual cmonkey2 runs of the ensemble have run, we need to
perform additional postprocessing.
Most of them are very compute-intensive task, which is why the tools actually
generate scripts to run on an SGE cluster.

Generate FIMO runs
-------------------

We have tested this with the fimo tool from the MEME suite 4.11.3. The command
line interface can changed between minor releases of MEME, so if you experience
problems, please verify that you have this version of MEME.

.. highlight:: none

::

  $ cd <target directory>
  $ egrin2-fimo --user <username> --organism <organism code> --genome <path pattern to genome files>

This generates the necessary scripts for running FIMO an SGE cluster. To submit
these jobs, on the head node, type

.. highlight:: none

::

  $ qsub qsub_fimo.sh

Coding Fractions
----------------

Immediately following the FIMO step is the Coding Fractions step. Coding fractions are
computed based on the FIMO results, so the entire FIMO step needs to be completed before
Coding Fractions can be run.

In order to generate scripts for running on the cluster, we need:

.. highlight:: none

::

  $ cd <target directory>
  $ egrin2-codingfracs --user <username> --organism <organism code> --features <RSAT features file>

and then on the SGE head node

.. highlight:: none

::

  $ qsub qsub_coding_fracs.sh



Tomtom
------

The TOMTOM step is somewhat independent of FIMO and Coding Fractions and could
theoretically be run in parallel if sufficient computing resources are available.
This command generates the control scripts for running TOMTOM.

.. highlight:: none

::

  $ cd <target directory>
  $ egrin2-tomtom --user <username> --prefix <output directory prefix>

then execute

.. highlight:: none

::

  $ qsub qsub_tomtom.sh

on the head node of your cluster.

Motif Clustering
----------------
