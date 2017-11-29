# vi: sw=4 ts=4 et:
"""datamatrix.py - data matrix/data frame manipulation for ensemble input data

This file is part of egrin2-tools. Please see README and LICENSE for
more information and licensing details.
"""
import scipy
import numpy as np
import pandas
import os
import random

# Python2/Python3 compatibility
try:
    xrange
except NameError:
    xrange = range


FILTER_THRESHOLD = 0.98
ROW_THRESHOLD = 0.17
COLUMN_THRESHOLD = 0.1


def nochange_filter(dataframe):
    """returns a new filtered DataMatrix containing only the columns and
    rows that have large enough measurements"""

    def nochange_filter_rows():
        """subfunction of nochange_filter to filter row-wise"""
        keep = []
        dmvalues = dataframe.values
        for row_index in xrange(dataframe.shape[0]):
            count = 0
            for col_index in xrange(dataframe.shape[1]):
                value = dmvalues[row_index, col_index]
                if np.isnan(value) or abs(value) <= ROW_THRESHOLD:
                    count += 1
            mean = float(count) / dataframe.shape[1]
            if mean < FILTER_THRESHOLD:
                keep.append(row_index)
        return keep

    def nochange_filter_columns():
        """subfunction of nochange_filter to filter column-wise"""
        keep = []
        dmvalues = dataframe.values
        for col_index in xrange(dataframe.shape[1]):
            count = 0
            for row_index in xrange(dataframe.shape[0]):
                value = dmvalues[row_index, col_index]
                if np.isnan(value) or abs(value) <= COLUMN_THRESHOLD:
                    count += 1
            mean = float(count) / dataframe.shape[0]
            if mean < FILTER_THRESHOLD:
                keep.append(col_index)
        return keep

    rows_to_keep = nochange_filter_rows()
    cols_to_keep = nochange_filter_columns()
    colnames = list(map(lambda col: dataframe.columns[col], cols_to_keep))
    rownames = list(map(lambda row: dataframe.index[row], rows_to_keep))
    numrows = len(rows_to_keep)
    numcols = len(cols_to_keep)

    rvalues = np.zeros((numrows, numcols))
    mvalues = dataframe.values
    for row_index in xrange(numrows):
        for col_index in xrange(numcols):
            value = mvalues[rows_to_keep[row_index], cols_to_keep[col_index]]
            rvalues[row_index, col_index] = value
    result = pandas.DataFrame(rvalues, rownames, colnames)
    return result


def row_filter(dataframe, fun):
    """generalize a matrix filter that is applying a function for each row"""
    num_rows = dataframe.shape[0]
    values = np.zeros(dataframe.shape)
    for row_index in xrange(num_rows):
        values[row_index] = fun(dataframe.values[row_index])

    return pandas.DataFrame(values, dataframe.index, dataframe.columns)


def r_stddev(values):
    """This is a standard deviation function, adjusted so it will
    return approximately the same value as R's sd() function would"""
    values = np.array(values)
    masked = values[np.isfinite(values)]
    num_values = len(masked)
    return round(np.std(masked) / np.sqrt(float(num_values - 1) /
                                          float(num_values)), 8)

def center_scale_filter(matrix):
    """center the values of each row around their median and scale
    by their standard deviation"""

    def center_scale(row):
        """centers the provided row around the median"""
        filtered = row[np.isfinite(row)]
        center = scipy.median(filtered)
        scale = r_stddev(filtered)
        nurow = [((value - center) / scale)
                 if not np.isnan(value) else value for value in row]
        return nurow

    return row_filter(matrix, center_scale)


def create_from_dataframe(dataframe, filters=[]):
    """creates and returns an initialized, filtered DataMatrix instance"""
    for matrix_filter in self.filters:
        data_matrix = matrix_filter(data_matrix)
        return data_matrix.sorted_by_row_name()

def prepare_ensemble_matrix(path, outdir, n, kmin):
    if os.path.exists(path):
        df = pandas.read_csv(path, index_col=0)
        for f in [nochange_filter, center_scale_filter]:
            df = f(df)
    return df

def split_matrix(df, outdir, n, kmin, kmax):
    """Split the input dataframe into n data frames  with the original
    number of rows and k columns. Write the resulting data frame to
    the specified output directory"""
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for i in xrange(n):
        k = random.randint(kmin, kmax)
        column_names = random.sample(list(df.columns), k)
        df_i = df[column_names]
        path = '%s/ratios-%03d.tsv' % (outdir, i + 1)
        df_i.to_csv(path)
