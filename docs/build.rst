Build
=====

Do you need an **EGRIN 2.0** ensemble ?
---------------------------------------

Good! You've come to the right place.

We've already developed ensembles for several commonly studied microbes, including E. coli and Mycobacterium tuberculosis. Consult our list of available ensembles to see if we've already developed a model for your organism.

If you don't see your organism listed, it's time to think about building an ensemble of your own. The tools described here will help you do that provided you have or can assemble a compendium of transcriptome profiles (e.g hundreds of microarrays). Larger datasets are generally better. If you need help, you can `contact us <https://www.systemsbiology.org/people/labs/baliga-lab/>`_ for guidance.

If you need to build an ensemble you should first familiarize yourself with the python-based `cMonkey2 <http://baliga-lab.github.io/cmonkey2/>`_ algorithm.

What is an **EGRIN 2.0** ensemble ?
------------------------------------

An EGRIN 2.0 ensemble is the result of running cMonkey2 (a non-deterministic, integrated biclustering algorithm) many times to infer highly-accurate, condition-specfic gene regulatory networks.

Why? Read our `paper <http://msb.embopress.org/content/10/7/740.long>`_.

What do the **BUILD** tools do ?
--------------------------------

These tools help build custom ensembles. In short, they automatically configure cMonkey2 runs and generate QSub scripts for submitting many ``cMonkey2`` jobs to a Sun Grid Engine (SGE) scheduler.

The tools generate cMonkey2 runs with slightly different parameterizations for training on subsets of the available experimental data. If you supply ``BLOCK`` definitions for the experimental conditions, including ``INCLUSION BLOCKS`` and ``EXCLUSION BLOCKS``, the subsets of experimental data can be chosen intelligently. More about this in the :doc:`build_tutorial`.

Requirements
------------
Obviously for these to work, you will need:

  * Access to a cluster running the SGE scheduler
  * `cMonkey2 <https://github.com/baliga-lab/cmonkey2>`_ installed on that cluster.
  * in addition, see system requirements listed in this documentation

Tutorial
--------
If you meet all of the requirements and are ready to build an ensemble from scratch, continue to our :doc:`build_tutorial`.
