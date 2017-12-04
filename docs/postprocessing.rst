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
  $ egrin2-fimo --user <username> --organism <organism code> --genome <path to genome file>

Coding Fractions
----------------

Tomtom
------

Motif Clustering
----------------
