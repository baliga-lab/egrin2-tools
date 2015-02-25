#!/usr/bin/env python

"""Tools for querying EGRIN2.0 MongoDB."""

__author__ = "Aaron Brooks"
__copyright__ = "Copyright 2014, cMonkey2"
__credits__ = ["Aaron Brooks", "David Reiss"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Aaron Brooks"
__email__ = "brooksan@uw.edu"
__status__ = "Development"

import random

from pymongo import MongoClient
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import itertools
from bson.code import Code
import matplotlib.pyplot as plt
import plotly.plotly as py
from plotly.graph_objs import *
import colorbrewer as cb
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

from assemble.resample import *
from query.egrin2_query import *

def plotExpression( data, plot_type = "boxplot", ipynb = False, zlim = None, sort = True, boxpoints = None ):
	
	if sort:
		c_order = data.mean(0).order().index.tolist()
		data = data.loc[:,c_order]

	def to_scatter( df ):
		x = df.index.values
		lines={}
		for key in df:
		    lines[key]={}
		    lines[key]["x"]=x
		    lines[key]["y"]=df[key].values
		    lines[key]["name"]=key

		    #Appending all lines
		lines_plotly=[lines[key] for key in df]
		return lines_plotly

	def to_box( df, boxpoints ):

		if boxpoints:
			boxpoints = "all"
		elif df.shape[1]<=50 and boxpoints is None:
			boxpoints = "all"
		elif boxpoints is None:
			boxpoints = False
		else:
			boxpoints = boxpoints

		boxes = []
		for x in df.columns.tolist():
			boxes.append(Box(
					name = x,
					y = df.loc[:,x].tolist(),
					boxpoints=boxpoints,
        			jitter=0.3,
        			pointpos=0
				))
		return boxes

	def to_heatmap( df, zlim ):

		if zlim is None:
			# use 1st and 99th percentile
			all_data = list( itertools.chain( *df.values.tolist() ) )
			if np.percentile( all_data, 1 ) < 0:
				zmin = -max( np.abs( np.percentile( all_data, [ 1, 99 ] ) ) )
			else:
				zmin = np.percentile( all_data, 1 )
			zmax = max( np.abs( np.percentile( all_data, [ 1, 99 ] ) ) )
		else:
			zmin = zlim[0]
			zmax = zlim[1]

		# do hierarchical clustering to make it look pretty
		D1 = squareform(pdist(df, metric='euclidean'))
		D2 = squareform(pdist(df.T, metric='euclidean'))
		Y = linkage(D1, method='complete')
		Z1 = dendrogram(Y, orientation='right')
		Y = linkage(D2, method='complete')
		Z2 = dendrogram(Y)
		idx1 = Z1['leaves']
		idx2 = Z2['leaves']
		D = df.iloc[idx1, :]
		D = D.iloc[:, idx2]
		to_r = Heatmap(
		        z = [ D.loc[i].tolist() for i in D.index.tolist() ],
		        y = D.index.tolist(),
		        x = D.columns.tolist(),
		        colorscale=[ [0.0, 'rgb(0,0,255)'], [0.5, 'rgb(0,0,0)'], [1.0, 'rgb(255,255,0)'] ],
		        zauto = False,
		        zmin = zmin,
		        zmax = zmax

		    )
		return to_r

	if plot_type == "line":
		to_plot = to_scatter(data.T)

		layout = Layout(
		title='Expression',
		xaxis=XAxis(
			title='Condition',
			ticks='',
			showticklabels=False
			),
		yaxis=YAxis(
			title='Expression Value',
			zeroline=True,
			)
		)

		fig = Figure(data=to_plot, layout=layout)

	elif plot_type == "heatmap":
		to_plot = to_heatmap(data, zlim)
		fig = Data( [ to_plot ] )
	elif plot_type == "boxplot":
		to_plot = to_box(data, boxpoints)
		fig = to_plot
	else:
		print "ERROR: Cannot recognize plot_type = %s" % plot_type
		return None


	if not ipynb:
		unique_url = py.plot( fig )
		return unique_url

	return fig