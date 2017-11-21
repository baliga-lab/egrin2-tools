Basic Mining Using ``agglom``
=============================

EGRIN2.0 MongoDB query using agglom
------------------------------------

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

``agglom`` function
~~~~~~~~~~~~~~~~~~~

The agglom function is the centerpiece for querying information in the ensemble.

Using this function you can co-associate any two bicluster features across the entire ensemble: so - biclusters, genes, conditions, and GREs with biclusters, genes, conditions, or GREs.

If you stack multiple agglom calls you can quickly build very complex queries.

The agglom function has a simple, universal format. Basically it needs an input (x), the type of that input (x_type), the desired output type (y_type), and some infomration about the location of the EGRIN 2.0 MongoDB.

The function will output a table of pvalues that reports the signifigance for each element in your output query type.

Let's take a look at the function itself, then we will explore serveral examples:

.. highlight:: none

::

   ?e2q.agglom


**Example 1: find co-regulated genes**

In this example, we will mine the ensemble for genes discovered in biclusters with the gene E. coli carbamoyl phosphate synthetase, carA (b0032)

.. highlight:: none

::

   carA_genes = e2q.agglom(db, x="b0032", x_type="genes", y_type="genes")
   carA_genes


===== ====== ========== ============= ============= ===============
gene  counts all_counts pval          qval_BH       qval_bonferroni
===== ====== ========== ============= ============= ===============
b1062 152    198        0.000000e+00  0.000000e+00  0.000000e+00
b0336 153    198        0.000000e+00  0.000000e+00  0.000000e+00
b0337 154    198        0.000000e+00  0.000000e+00  0.000000e+00
b4244 167    198        0.000000e+00  0.000000e+00  0.000000e+00
b4245 169    198        0.000000e+00  0.000000e+00  0.000000e+00
b0033 173    198        0.000000e+00  0.000000e+00  0.000000e+00
b0032 198    198        0.000000e+00  0.000000e+00  0.000000e+00
b2499 142    198        3.783051e-300 4.073085e-298 3.665776e-297
b2476 142    198        3.783051e-300 4.073085e-298 3.665776e-297
b2313 141    198        5.664712e-297 4.574255e-295 5.489106e-294
b2312 141    198        5.664712e-297 4.574255e-295 5.489106e-294
...   ...    ...        ...           ...           ...
===== ====== ========== ============= ============= ===============

259 rows × 5 columns

There are several things to note in this query.

First the arguments:

  * x specfies the query. This can be gene(s), condition(s), GRE(s), or bicluster(s). x can be a single entity or a list of entitites of the same type.
  * x_type indicates the type of x. This can include gene, condition, gres, and biclusters. Basically: "what is x?" The parameter typing is pretty flexible, so - for example - rows can be used instead of genes.
  * y_type is the type of. Again, genes, conditions, gres, or biclusters.
  * host specifies where the MongoDB database is running. In this case it is running on a machine called baligadev. If you are hosting the database locally this would be localhost
  * db is the name of the database you want to perform the query in. Typically databases are specified by the three letter organism code (e.g., eco) followed by _db. A list of maintained databases is available here.

Also note the output:

This is a table of all genes that are frequently co-discovered in biclusters with b0032 (carA). They are sorted by their relative "enrichment" in the sample x (where lower pval = greater enrichment). This is quantified by the hypergeometric distribution. To account for multiple testing, two corrections are included: Benjamini-Hochberg correction (qval_BH; less stringent) and Bonferrnoni correction (qval_bonferroni; more stringent).

**Example 2: what if gene names have different formats ?**

We can easily change the format of x, which makes it easy to specify genes in any format, e.g. NCBI annotation or common name. Let's try to include a list of genes with two different naming formats. This time we will only return the first few elements using .head()

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type= "genes").head()

===== ====== ========== ==== ======= ===============
gene  counts all_counts pval qval_BH qval_bonferroni
===== ====== ========== ==== ======= ===============
b4244 179    198        0    0       0
b0033 180    198        0    0       0
b4245 185    198        0    0       0
b0031 198    198        0    0       0
b4246 198    198        0    0       0
===== ====== ========== ==== ======= ===============

5 rows × 5 columns

Here, two of the genes are Blattner Ids (b****), while the other is a common name. This is no problem for the query.

The lookup table is pretty extensive for genes (it comes from MicrobesOnline) so anything from common name to GI number can be matched. Cool, huh?

**Example 3: alternative logic**

You can specify the type of logic to use in grouping your query.

For example, do you want do find biclusters where all genes x in your query are present (and logic) or biclusters where at least one of the genes is present (or logic). You can change the logical grouping of your query by specifying the logic parameter.

This parameter can take one of three values: and, or, or nor. By default the command uses or logic.

Here we change the logic employed in the previous query to "and", meaning that a bicluster must contain all of the query genes. Notice how the values change.

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type="genes", logic="and").head()

===== ====== ========== ============ ============ ===============
gene  counts all_counts pval         qval_BH      qval_bonferroni
===== ====== ========== ============ ============ ===============
b4246 25     198        0.000000e+00 0.000000e+00 0.000000e+00
b0032 25     198        0.000000e+00 0.000000e+00 0.000000e+00
b0031 25     198        0.000000e+00 0.000000e+00 0.000000e+00
b4006 24     198        5.510260e-57 3.981163e-55 1.592465e-54
b1062 23     198        2.604624e-53 1.075338e-51 7.527363e-51
===== ====== ========== ============ ============ ===============

5 rows × 5 columns

Logical operations can provide a substantial amount of power to your queries. For example, imagine that you not only want to know when b0032, pyrL, and b0031 are co-regulated but also when they are not. This can be done using the logical operator nor.

**Example 4: other associations**

So far we have mined gene-gene co-associations, but there are other associations in the ensemble that are also available, including conditions, and GREs, and biclusters.

Each of these can be accessed by changing the x_type and y_type parameter.

**Conditions**

For example to find conditions where these genes are co-regulated, we change y_type to conditions:

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type="conditions", logic="or").head()

======================= ====== ========== ============ ============ ===============
conditions              counts all_counts pval         qval_BH      qval_bonferroni
======================= ====== ========== ============ ============ ===============
rb_del_rpoS_biofilm     240    12921      2.849586e-13 7.228851e-11 1.327907e-10
biofilm_24hr_del_yceP   239    12879      3.889790e-13 7.228851e-11 1.812642e-10
ik_H2_T4                232    12380      4.653767e-13 7.228851e-11 2.168655e-10
MG1063_uninduced_t180   232    12543      2.190004e-12 2.551355e-10 1.020542e-09
rb_del_rpoS_exponential 226    12198      5.366065e-12 5.001173e-10 2.500586e-09
======================= ====== ========== ============ ============ ===============

5 rows × 5 columns

**GREs**

To discover GREs located upstream of these genes we change y_type to gre:

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type="gre", logic="or").head()

==== ====== ========== ============= ============= ===============
GREs counts all_counts pval          qval_BH       qval_bonferroni
==== ====== ========== ============= ============= ===============
11   187    231        9.726498e-267 5.252309e-265 5.252309e-265
81   42     48         2.936386e-64  7.928242e-63  1.585648e-62
12   41     211        6.248311e-26  1.124696e-24  3.374088e-24
309  14     19         2.195152e-21  2.606271e-20  1.185382e-19
197  16     27         2.413214e-21  2.606271e-20  1.303136e-19
==== ====== ========== ============= ============= ===============

5 rows × 5 columns

Note that GRE names in ensemble are integer values, thus GRE #1 is simply called 1

**Biclusters**

Finally: to return biclusters we change y_type to bicluster.

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type="bicluster", logic="or").head()

   WARNING:root:Will return bicluster _id. The results might surprise you...

========= ========================
Bicluster _id
========= ========================
0         54ee42396a208812a205295a
1         54ee42396a208812a2052982
2         54ee42826a208812a2052b61
3         54ee42836a208812a2052ce0
4         54ee430d6a208812a2053120
========= ========================

5 rows × 1 columns

Please note that this returns MongoDB _ids for the biclusters, which are not by themselves useful or informative. Rather, these _ids can be used in other queries (e.g. fimoFinder) or to query MongoDB directly.
