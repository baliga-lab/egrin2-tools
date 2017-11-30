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
