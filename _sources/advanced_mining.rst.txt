Advanced Mining
===============

In a nutshell
-------------

Information hosted in the EGRIN 2.0 MongoDB databases can be used directly to extend EGRIN 2.0 queries and analysis.

Set-up
------

As usual, make sure ``./egrin-tools/` folder is in your ``$PYTHONPATH``.
You can do this on Mac/Linux by adding the path to you Bash Shell Startup Files, e.g. ``~/.bashrc`` or ``~/.bash_profile``

for example in ``~/.bash_profile`` add the following line:

.. highlight:: none

::

   export PYTHONPATH=$PYTHONPATH:path/to/egrin2-tools/

Load required modules
---------------------

.. highlight:: none

::

   import query.egrin2_query as e2q
   import pymongo
   import pandas as pd

   client = pymongo.MongoClient(host='localhost')
   db = client['eco_db']


There are several dependencies that need to be satisfied, including:

  * ``pymongo``
  * ``numpy``
  * ``pandas``
  * ``joblib``
  * ``scipy``
  * ``statsmodels``
  * ``itertools``

Run agglom function
-------------------

**Example 1: Find genes co-regulated in biclusters with the gene carA (b0032)**

.. highlight:: none

::

   carA_genes = e2q.agglom(db, x="b0032", x_type="genes", y_type="genes")
   carA_genes

===== ====== ========== ============= ============= ===============
gene  counts all_counts pval          qval_BH       qval_bonferroni
===== ====== ========== ============= ============= ===============
b0032 198    198        0.000000e+00  0.000000e+00  0.000000e+00
b0033 173    198        0.000000e+00  0.000000e+00  0.000000e+00
b4245 169    198        0.000000e+00  0.000000e+00  0.000000e+00
b4244 167    198        0.000000e+00  0.000000e+00  0.000000e+00
b0337 154    198        0.000000e+00  0.000000e+00  0.000000e+00
b0336 153    198        0.000000e+00  0.000000e+00  0.000000e+00
b1062 152    198        0.000000e+00  0.000000e+00  0.000000e+00
b2476 142    198        3.783051e-300 4.073085e-298 3.665776e-297
b2499 142    198        3.783051e-300 4.073085e-298 3.665776e-297
b2312 141    198        5.664712e-297 4.574255e-295 5.489106e-294
b2313 141    198        5.664712e-297 4.574255e-295 5.489106e-294
b2557 141    198        5.664712e-297 4.574255e-295 5.489106e-294
b0522 138    198        1.477589e-287 1.101372e-285 1.431783e-284
b1849 136    198        2.280977e-281 1.578762e-279 2.210267e-278
b4006 132    198        3.382871e-269 2.185335e-267 3.278002e-266
b0523 129    198        3.043903e-260 1.843464e-258 2.949542e-257
b3654 128    198        2.728256e-257 1.555106e-255 2.643680e-254
b4246 126    198        1.964952e-251 1.057799e-249 1.904038e-248
b2500 122    198        6.646122e-240 3.389522e-238 6.440093e-237
b2497 112    198        4.061863e-212 1.967973e-210 3.935945e-209
b4064 102    198        9.823972e-186 4.533061e-184 9.519429e-183
b0945 96     198        1.511792e-170 6.658758e-169 1.464927e-167
b4005 89     198        2.020583e-153 8.512805e-152 1.957945e-150
b3941 58     198        3.940658e-85  1.591041e-83  3.818498e-82
...   ...    ...        ...           ...           ...
===== ====== ========== ============= ============= ===============

259 rows × 5 columns

There are several things to note in this query.

First the arguments:

  * ``x`` specfies the query. This can be gene(s), condition(s), GRE(s), or bicluster(s). x can be a single entity or a list of entitites of the same type.
  * ``x_type`` indicates the type of x. This can include gene, condition, gres, and biclusters. Basically: what is x? The parameter typing is pretty flexible, so - for example - rows can be used instead of genes.
  * ``y_type`` is the type of. Again, genes, conditions, gres, or biclusters.
  * ``host`` specifies where the MongoDB database is running. In this case it is running on a machine called primordial (my box). If you are hosting the database locally this would be localhost
  * ``db`` is the name of the database you want to perform the query in. Typically databases are specified by the three letter organism code (e.g., eco) followed by _db.

Also note the output:

This is a table of all genes that are frequently co-discovered in biclusters with b0032 (carA). They are sorted by their relative "enrichment" in the sample x (where lower pval = greater enrichment). This is quantified by the hypergeometric distribution. To account for multiple testing, two corrections are included: Benjamini-Hochber correction (qval_BH; less stringent) and Bonferrnoni correction (qval_bonferroni; more stringent).

We can easily change the format of x, which makes it easy to specify genes in any format. Let's try to include a list of genes with two different naming formats. This time we will only return the first few elements using .head()


.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type="genes").head()

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

The lookup table is pretty extensive for genes (it comes from MicrobesOnline) so anything from common name to GI number can be matched. Cool, huh?

You can also specify the type of logic to use in grouping your query. For example, do you want do find biclusters where all genes x in your query are present (and logic) or biclusters where at least one of the genes is present (or logic). You can change the logical grouping of your query by specifying the logic parameter. This parameter can take one of three values: ``and``, ``o``r, or ``nor``. By default the command uses and logic.

Here we change the logic employed in the previous query. Notice how the values change.

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type="genes", logic="or").head()

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

Logical operations can provide a substantial amount of power to your queries. For example, imagine that you not only want to know when b0032, pyrL, and b0031 are co-regulated but also when they are not. This can be done using the logical operator nor.

First we query for conditions where they are co-regulated:

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type="conditions", logic="or").head()


======================= ============ ============ ============ ============ ===============
condition               counts       all_counts   pval         qval_BH      qval_bonferroni
======================= ============ ============ ============ ============ ===============
rb_del_rpoS_biofilm     240          12921        2.849586e-13 7.228851e-11 1.327907e-10
biofilm_24hr_del_yceP   239          12879        3.889790e-13 7.228851e-11 1.812642e-10
ik_H2_T4                232          12380        4.653767e-13 7.228851e-11 2.168655e-10
MG1063_uninduced_t180   232          12543        2.190004e-12 2.551355e-10 1.020542e-09
rb_del_rpoS_exponential 226          12198        5.366065e-12 5.001173e-10 2.500586e-09
======================= ============ ============ ============ ============ ===============

5 rows × 5 columns

Now we can compare these conditions to conditions that are exluded from biclusters where at least one of the genes is present:

.. highlight:: none

::

   e2q.agglom(db, x=["b0032","pyrL","b0031"], x_type="genes", y_type= "conditions", log


========================================= ====== ========== ============ ======== ===============
condition                                 counts all_counts pval         qval_BH  qval_bonferroni
========================================= ====== ========== ============ ======== ===============
har_S0_R_noIPTG                           12933  13045      2.442307e-09 0.000001 0.000001
MOPS_K_dps_stationary2                    10149  10230      5.916970e-09 0.000001 0.000003
BW25113_K_tnaA_30C_biofilm_indole_control 11308  11403      9.537266e-09 0.000001 0.000004
MG1655_kanamycin_t60                      11240  11337      5.362424e-08 0.000006 0.000025
dilution_t3                               12110  12219      1.446307e-07 0.000013 0.000067
========================================= ====== ========== ============ ======== ===============

5 rows × 5 columns

Understanding condition-specific co-regulation of gene modules is what EGRIN 2.0 is all about. Because MongoDB allows for flexible document structures, we've included options for rich experimental annotation. Providing details about what is going on in an experiment (metainformation) provides an additional layer of inquiry (and hopefully insight). In the next few examples, I demonstrate how you might use the experimental annotations in interesting ways.

First, let's get an idea for what these experimental annotations might look like. Let's use the condition called rb_del_rpoS_biofilm that came up as a condition where b0032, pyrL, and b0031 are co-regulated.

To do this we need to connect to the database.

.. highlight:: none

::

   client = pymongo.MongoClient(host="localhost", port=27017 )

Now we can query the col_info collection in our database. Here I am just returning the field called additional_info, which is where we store all of the optional annotations. These annotations come directly from the M3D database. The full databse schema is available for your reference on the wiki.

.. highlight:: none

::

   client["eco_db"].col_info.find_one({"egrin2_col_name": "rb_del_rpoS_biofilm"}, {"_id": 0, "additional_info": 1})
   {u'additional_info': [{u'name': u'aeration',
      u'units': nan,
      u'value': u'assumed_anaerobic'},
     {u'name': u'ammonium_chloride', u'units': u'mM', u'value': u'9.52'},
     {u'name': u'ammonium_molybdate', u'units': u'mM', u'value': u'0.00000291'},
     {u'name': u'boric_acid', u'units': u'mM', u'value': u'0.000401'},
     {u'name': u'calcium_chloride', u'units': u'mM', u'value': u'0.0005'},
     {u'name': u'cobalt_chloride', u'units': u'mM', u'value': u'0.0000303'},
     {u'name': u'culture_temperature', u'units': u'Celsius', u'value': u'37'},
     {u'name': u'culture_type', u'units': nan, u'value': u'fed_batch'},
     {u'name': u'culture_vessel', u'units': nan, u'value': u'flow cell'},
     {u'name': u'cupric_sulfate', u'units': u'mM', u'value': u'0.00000961'},
     {u'name': u'experimenter', u'units': nan, u'value': u'Ito A'},
     {u'name': u'ferrous_sulfate', u'units': u'mM', u'value': u'0.01'},
     {u'name': u'glucose', u'units': u'mM', u'value': u'11.1'},
     {u'name': u'growth_phase', u'units': nan, u'value': u'biofilm'},
     {u'name': u'magnesium_chloride', u'units': u'mM', u'value': u'0.525'},
     {u'name': u'manganese(II)_chloride',
      u'units': u'mM',
      u'value': u'0.0000808'},
     {u'name': u'MOPS', u'units': u'mM', u'value': u'40'},
     {u'name': u'note',
      u'units': nan,
      u'value': u'channel dimensions: 1x4x40 mm, medium pumped thru at 0.25 mL/min, samples taken at 72 hrs'},
     {u'name': u'perturbation', u'units': nan, u'value': u'knockout'},
     {u'name': u'perturbation_gene', u'units': nan, u'value': u'rpoS'},
     {u'name': u'potassium_phosphate_monobasic',
      u'units': u'mM',
      u'value': u'1.32'},
     {u'name': u'potassium_sulfate', u'units': u'mM', u'value': u'0.276'},
     {u'name': u'RNA_prep_type', u'units': nan, u'value': u'RNeasy'},
     {u'name': u'RNA_stop_solution', u'units': nan, u'value': u'RNAprotect'},
     {u'name': u'sodium_chloride', u'units': u'mM', u'value': u'50'},
     {u'name': u'strain', u'units': nan, u'value': u'MG1655 del rpoS'},
     {u'name': u'structured_metadata', u'units': nan, u'value': u'complete'},
     {u'name': u'thiamine_HCl', u'units': u'mM', u'value': u'0.2964'},
     {u'name': u'tricine', u'units': u'mM', u'value': u'4'},
     {u'name': u'zinc_sulfate', u'units': u'mM', u'value': u'0.00000974'}]}

Notice that each optional annotation has name, units, and value fields. We can take advatage of this information and the nice query features of MongoDB to compose more powerful queries.

For example, you might have noticed in our original query that several conditions with biofilm in their titles came up. Maybe our genes play some role in biofilm-related processes? Let's see if we can return all genes that are co-regulated in conditions where growth_phase is biofim.

To do that we need to retrieve all of the conditions where growth_phase is biofilm.

.. highlight:: none

::

   biofilm_conds = pd.DataFrame(list(client["eco_db"].col_info.find({"additional_info.name": "growth_phase", "additional_info.value": "biofilm"}, {"_id": 0, "egrin2_col_name": 1})))

   print biofilm_conds.head()

     egrin2_col_name
   0 biofilm_K_yceP
   1 biofilm_K_yceP_indole
   2 biofilm_wt_glucose
   3 biofilm_K_trpE
   4 biofilm_K_tnaA

   [5 rows x 1 columns]

   biofilm_conds.count()

   egrin2_col_name    46
   dtype: int64

There are 46 conditions annotated as biofilm. Just for fun, let's say we wanted restrict these conditions to those where the strain used was *E. coli K-12 MG1655.*

.. highlight:: none

::

   biofilm_conds_MG1655 = pd.DataFrame(list(client["eco_db"].col_info.find({"$and": [{"additional_info.name": "strain", "additional_info.value": {"$regex": "MG1655" } }, { "additional_info.name": "growth_phase", "additional_info.value": "biofilm"}]}, {"_id":0, "egrin2_col_name": 1})))

   biofilm_conds_MG1655.head()

== ==============================
#  egrin2_col_name
== ==============================
0  MG1655_wt_24hr_biofilm
1  MG1655_wt_R1drd19_24hr_biofilm
2  rb_del_rpoS_biofilm
3  rb_wt_biofilm
4  bform_biofilm_attachment
== ==============================

.. highlight:: none

::

   biofilm_conds_MG1655.count()
   egrin2_col_name    11
   dtype: int64

You can see that these 11 are a subset of the 46 originally returned. This query highlights some features of MongoDB queries. Here we've used a regular expression to find documents where the strain annotation contained the values "MG1655". You can also find documents by value comparisons (e.g. equalities/inequalities). The document searching capabilities of MongoDB are pretty extensive. You can learn more about query operators in MongoDB here.

Now we can find the genes that are co-regulated in the 11 experiments where the growth phase was biofilm. Here I will use the and logical operator, which would be the most stringent criteria (i.e. all conditions must be present in the bicluster).

.. highlight:: none

::

   biofilm_genes = e2q.agglom(db, x=biofilm_conds_MG1655.egrin2_col_name.tolist(), x_type="conditions", y_type="genes", logic="and")

   print biofilm_genes.head()

===== ====== ========== ======== ======== ===============
gene  counts all_counts pval     qval_BH  qval_bonferroni
===== ====== ========== ======== ======== ===============
b3963 4      198        0.000005 0.001972 0.003943
b2796 4      198        0.000005 0.001972 0.003943
b3774 3      198        0.000116 0.002955 0.088654
b4443 3      198        0.000116 0.002955 0.088654
b3357 3      198        0.000116 0.002955 0.088654
===== ====== ========== ======== ======== ===============

[5 rows x 5 columns]

We can restrict this list to genes that surpass a multiple testing correction at a qval threshold of 0.05.

.. highlight:: none

::

   biofilm_genes_sig = biofilm_genes.qval_BH[biofilm_genes.qval_BH < 0.05]
   biofilm_genes_sig.shape[0]
   761

There are 761 genes that meet our criteria. Now we can see if our genes are in this list.

.. highlight:: none

::

   "b0032" in biofilm_genes.qval_BH[biofilm_genes.qval_BH < 0.05]
   True

   e2q.row2id_batch(db, ["pyrL"], return_field = "egrin2_row_name")[0] in biofilm_genes.qval_BH[biofilm_genes.qval_BH < 0.05]
   False

   "b0031" in biofilm_genes.qval_BH[biofilm_genes.qval_BH < 0.05]
   False

So 1 out of 3 of our genes occur frequently in biolfim annotated conditions. Not very strong support for our hypothesis.

In this example, however, we also see an example of name translation using the row2id_batch() function. Here we need to translate the common name pyrL to the name used in EGRIN 2.0, b4246. We did that as follows:

.. highlight:: none

::

   print e2q.row2id_batch(db, ["pyrL"], return_field="egrin2_row_name")[0]
   b4246
