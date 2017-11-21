Basic Mining Using ``coremFinder``
==================================

EGRIN2.0 MongoDB query using coremFinder
----------------------------------------

In a nutshell
~~~~~~~~~~~~~

Corems or condition-specific co-regulated modules are sets of genes that are tightly co-expressed in a condition-specific manner and whose expression is often controlled by common transcriptional regulators.

Expert note: A different perspective on corems is that they are highly-refined, reproducibly-detected, biclusters that violate some constraints imposed on cMonkey-detected biclusters.

coremFinder mines information about corems that are detected by EGRIN 2.0.

Typically, about half of the genes in the genome may be discovered in corems. For this subset of genes, the coremFinder function can return information about the corem, including:

  * gene composition
  * condition-specific activity
  * edges contained in corem

This function can be combined with the agglom function to drive analysis of EGRIN 2.0 ensembles beyond genes contained in corems.

Set-up
~~~~~~

Make sure ``./egrin-tools/`` folder is in your ``$PYTHONPATH``

You can do this on Mac/Linux by adding the path to you Bash Shell Startup Files, e.g. ``~/.bashrc`` or ``~/.bash_profile``

for example in ``~/.bash_profile`` add the following line:

.. highlight:: none

::

   export PYTHONPATH=$PYTHONPATH:path/to/egrin2-tools/

Load required modules
~~~~~~~~~~~~~~~~~~~~~

.. highlight:: none

::

   import query.egrin2_query as e2q
   import pymongo

   client = pymongo.MongoClient(host='baligadev')
   db = client['eco_db']

There are several dependencies that need to be satisfied, including:

  * pymongo
  * numpy
  * pandas
  * joblib
  * scipy
  * statsmodels
  * itertools

find_corem_info() function
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``find_corem_info()`` function is very similar to the agglom() function. Basically, it co-associates information about corems, where the infomration supplied and retrieved is modulated by defining the arguments: ``x``, ``x_type``, and ``y_type``.

The function returns the requested information about a corem.

You can find out more about this function and its parameters by issuing the following commmand:

.. highlight:: none

::

   ?e2q.find_corem_info

**Example 1: Find corem genes**

The most straightforward way to use coremFinder is to find genes contained in a corem.

For example, to find all genes in E. coli corem #1 we would type:

.. highlight:: none

::

   corem_1 = e2q.find_corem_info(db, x=1, x_type="corem", y_type="genes")
   corem_1

===== =======
corem genes
===== =======
0     b3317
1     b3320
2     b3319
3     b3313
4     b3315
5     b3318
6     b3314
7     b3321
8     b3316
===== =======

9 rows × 1 columns

There are several things to note in this query.

First the arguments:

  * ``x`` specfies the query. This can be gene(s), condition(s), GRE(s), or edge(s). x can be a single entity or a list of entitites of the same type.
  * ``x_type`` indicates the type of x. This can include gene, condition, gres, and edges. Basically: "what is x?" The parameter typing is pretty flexible, so - for example - rows can be used instead of genes.
  * ``y_type`` is the type of. Again, genes, conditions, gres, or edges.
  * ``host`` specifies where the MongoDB database is running. In this case it is running on a machine called baligadev. If you are hosting the database locally this would be localhost
  * ``db`` is the name of the database you want to perform the query in. Typically databases are specified by the three letter organism code (e.g., eco) followed by _db. A list of maintained databases is available here.

Also notice that corems (like GREs) are named as integer values.

It should also be noted that corems are ordered by their weighted-density. Thus, corem #1 is the most densly connected corem in the network. Basically, this means that each gene in the corem is co-discovered frequently in biclusters with every other gene in that corem (strongly connected subnetwork).

Here we see that if we translate the names of these genes, we find that they are part of a ribosomal operon, which makes sense in light of the fact that ribosomal genes are tightly co-expressed.

.. highlight:: none

::

   e2q.row2id_batch(db, corem_1.genes.tolist(), return_field="name")

   [u'rplB',
    u'rplC',
    u'rplD',
    u'rplP',
    u'rplV',
    u'rplW',
    u'rpsC',
    u'rpsJ',
    u'rpsS']

**Example 2: Find corems for a specific gene**

More commonly, you want to know the corems to which a particular gene belongs.

This can be accomplished by changing x, x_type, and y_type, as follows:

.. highlight:: none

::

   carA_corems = e2q.find_corem_info(db, x="carA", x_type="gene", y_type="corems")
   carA_corems

= ======
# corems
= ======
0 107
1 471
2 835
3 847
= ======

4 rows × 1 columns

We can see from this query that carA belongs to four corems. We could retrieve the genes in these corems like in Example 1:

.. highlight:: none

::

   e2q.find_corem_info(db, x=carA_corems.corems.tolist(), x_type="corems", y_type="genes")

== =====
#  genes
== =====
0  b0002
1  b0003
2  b0004
3  b0032
4  b0033
5  b0197
6  b0198
7  b0273
8  b0287
9  b0336
10 b0337
.. ...
== =====

75 rows × 1 columns

**Example 3: Logical operations**

Similar to the agglom function we can implement logical operations. For example, if we wanted to know the genes that belonged to all of the corems in which carA is a member we would simply set logic = "and"

.. highlight:: none

::

   e2q.find_corem_info(db, x=carA_corems.corems.tolist(), x_type="corems", y_type="genes", logic="and")

== =====
#  genes
== =====
0  b0032
1  b0033
2  b0197
3  b0198
4  b0287
5  b2500
6  b2600
7  b2601
8  b3769
9  b3770
10 b3771
11 b3772
12 b3956
13 b4005
14 b4006
15 b4064
16 b4246
17 b4488
== =====

18 rows × 1 columns

Notice that only 17 out of the 75 genes in these four corems are present in every one of the four corems

**Example 4: Corem discovery based on experimental conditions**

Similar to the gene example above, corems can be discovered based on the experimental conditions in which the genes in a corem are co-regulated as well. For example, the experimental conditions associated with corem #1 can be discovered by changing the y_type supplied to one of the previous commands. Since there are many conditions, we will only display the first 10.

.. highlight:: none

::


   corem_1_conditions = e2q.find_corem_info(db, x=1, x_type="corem", y_type="conditions")
   corem_1_conditions[0:10]

== =====================
#  conditions
== =====================
0  str_ctrl_0m
1  str_str_K_relA_M9
2  str_ctrl_K_relA_M9
3  str_ctrl_M9
4  W3110_wt_luxS_glucose
5  W3110_K_luxS_glucose
6  suspension_24hr
7  suspension_15hr
8  biofilm_15hr
9  str_str_LV_20m
== =====================

10 rows × 1 columns

Likewise, we could retrieve all of the other corems that are also "active" in the first 10 conditions annotated to corem #1 by:

.. highlight:: none

::

   e2q.find_corem_info(db, x=corem_1_conditions.conditions[0:10].tolist(), x_type="conditions", y_type="corems", logic="and")

== ======
#  corems
== ======
0  1
1  3
2  45
3  46
4  55
5  74
6  76
7  104
8  111
9  114
10 117
.. ...
== ======

56 rows × 1 columns

Thus there are 55 corems that are co-regulated in all of the (first) ten conditions in which the genes in corem #1 are also co-regulated.

**Example 4: Edges**

Technically, corems are "link-communities", meaning that they are sets of edges, where the edge is a co-regulatory assocaition between two genes (nodes). This is why a single gene (node) can belong to multiple corems (link-communities).

To retrieve that actual edges that define a corem, set ``y_type`` to "``edges``":

.. highlight:: none

::

   e2q.find_corem_info(db, x=1, x_type="corem", y_type="edges")

== ===========
#  edges
== ===========
0  b3313-b3314
1  b3313-b3315
2  b3313-b3316
3  b3313-b3317
4  b3313-b3318
5  b3313-b3319
6  b3313-b3320
7  b3313-b3321
8  b3314-b3315
9  b3314-b3316
10 b3314-b3317
11 b3314-b3318
12 b3314-b3319
13 b3314-b3320
14 b3314-b3321
15 b3316-b3315
16 b3317-b3315
17 b3317-b3316
18 b3318-b3315
19 b3318-b3316
20 b3318-b3317
21 b3318-b3320
22 b3319-b3315
23 b3319-b3316
24 b3319-b3317
25 b3319-b3318
26 b3319-b3320
27 b3320-b3315
28 b3320-b3316
29 b3320-b3317
30 b3321-b3315
31 b3321-b3316
32 b3321-b3317
33 b3321-b3318
34 b3321-b3319
35 b3321-b3320
== ===========

36 rows × 1 columns
