Assemble
========

In a nutshell
-------------

The ``ASSEMBLE`` functions post-process individual cMonkey2 runs. These steps include:

  * GRE discovery by ``TOMTOM``
  * Genome-wide GRE scanning by ``FIMO``
  * Aggregation of all ensemble info into MongoDB
  * Corem detection

These tools return a MongoDB database that can be queried by the QUERY functions

Requirements
------------

For these scripts to work, you will need:

  * Access to a cluster running the SGE scheduler
  * A completed cMonkey2 ensemble (see BUILD documentation)

in addition to the system requirements listed in this documentation

Tutorial
--------

If you meet all of the requirements and are ready to ``ASSEMBLE`` an ensemble, continue to our :doc:`assemble_tutorial`.

MongoDB Database Schema
-----------------------

The MongoDB Database Schema is available here
