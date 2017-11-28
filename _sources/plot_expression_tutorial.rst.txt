TUTORIAL: Visualizing Basic `EGRIN 2.0` queries
===============================================

*In this tutorial, we will visualize gene expression resulting from basic EGRIN 2.0 queries*

Preliminaries
-------------

As described in the System Requirements, there are several dependencies that need to be satisfied to complete this tutorial, including:

  * ``pymongo``
  * ``numpy``
  * ``pandas``
  * ``joblib``
  * ``scipy``
  * ``statsmodels``
  * ``itertools``

To query the ensemble, we must first load all of the query functions.

.. highlight:: none

::

   import query.egrin2_query as e2q
   import query.egrin2_plot as e2p
   import pymongo


We will also define the host and the database that we'd like to use. Host is the name of the machine hosting the EGRIN 2.0 MongoDB while db is the organism-specific EGRIN 2.0 database name to query.

.. highlight:: none

::

   client = pymongo.MongoClient(host='baligadev')
   db = client['eco_db']

Basic queries
-------------

STEP 1: Find genes in a corem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we will retrieve genes and experiments in which these genes are co-regulated from a specific corem.

This type of information can be retrieved using the coremFinder function. To call this function we need to specify:

  * ``x``: our query
  * ``x_type``: our query type. This could be corems, genes, conditions, GREs, or specific-coregulatory edges. In this case we will use corem.
  * ``y_type``: our target type. This is the type of information we would like to retrieve. The type can be any type described by x_type

In addition we include the host and db variables defined above.

.. highlight:: none

::

   corem = 1
   corem_genes = e2q.find_corem_info(db, corem, x_type="corem_id", y_type="genes")
   corem_genes.sort()
   print "\nThere are %d genes in corem %d, including:\n" % (len(corem_genes), corem)
   for i in corem_genes:
       print i, " is also called ", e2q.row2id_batch(db, [i], return_field="name" )[0]

(TODO)


STEP 2: Find experiments where these genes are co-expressed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This query will also use the coremFinder function. To find conditions associated with our corem rather than genes, we simply change the y_type argument.

.. highlight:: none

::

   corem_conditions = e2q.find_corem_info(db, corem, x_type="corem", y_type="conditions")
   corem_conditions.sort()
   print "\nThere are %d conditions in which these genes are co-regulated, including:\n" % len(corem_conditions)
   for i in corem_conditions[0:10]:
       print i

   There are 418 conditions in which these genes are co-regulated, including:

   conditions


STEP 3: Retrieve gene expression from the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To retrieve gene expression values for these genes and experiments we use the expressionFinder function. To call this function we specific the rows (genes) and columns (conditions), as well as the host and db as before.

.. highlight:: none

::

   gene_expression = e2q.find_gene_expression(db, rows=corem_genes, cols=corem_conditions)

(TODO)


STEP 4. Plot expression values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have several options for plotting these gene expression values. We could plot the expression values as lines, or in a heatmap, or even as a boxplot for all genes in each condition.

Each of these visualizations is available by calling a single function, plotExpression. To call this function we must provide:

  * ``data``: this is the gene expression values, a Pandas data frame
  * ``plot_type``: this is the type of plot to draw. Can be ``boxplot`` (default),line, or ``heatmap``
  * ``ipynb``: logical indicating whether the plot will be drawn in an iPython notebook
  * ``sort``: optionally sort the data (default: ``FALSE``)

It's important to note that this function only requires a Pandas data frame, meaning you can use it to plot any kind of data (e.g. loaded from a text file).

Additionally, if we are producing the plot in an iPython notebook, we need to set the argument ipynb = TRUE and call an additional function py.iplot on the value returned from plotExpression.

Below are three examples, calling the plotExpression function with three different values for the plot_type argument: line, heatmap, and boxplot.

.. highlight:: none

::

   line_plot = e2p.plotExpression(data=gene_expression, plot_type = "line", ipynb=True, sort=False)
   py.iplot(line_plot)

(TODO)

.. highlight:: none

::

   heatmap = plotExpression( data = gene_expression, plot_type = "heatmap", ipynb = True )
   py.iplot( heatmap )

(TODO)

.. highlight:: none

::

   boxplot = plotExpression( gene_expression, plot_type = "boxplot", ipynb = True, sort = True )
   py.iplot( boxplot )

(TODO)

Putting it all together
~~~~~~~~~~~~~~~~~~~~~~~

Here is the code. You can copy this to your own notebook or download a blank slate notebook here.

.. highlight:: none

::

   # PLOT GENE EXPRESSION

   # prelims
   from query.egrin2_query import *
   host = ""
   db = ""

   corem = 1

   # find corem genes
   corem_genes =  coremFinder(x = corem,x_type = "corem", y_type="genes",host=host,db=db)
   # find corem conditions
   corem_conditions =  coremFinder(x = corem,x_type = "corem", y_type="conditions", host=host, db=db)
   # get gene expression
   gene_expression = expressionFinder(rows=corem_genes,cols=corem_conditions,host=host,db=db)
   # plot
   plot = plotExpression( data = gene_expression, plot_type = "line", ipynb = True, sort = False )
   py.iplot( plot )
