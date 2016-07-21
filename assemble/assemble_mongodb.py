import pandas as pd


class MongoDB:
    def __init__(self, dbclient):
        self.dbclient = dbclient

    def close(self):
        pass

    def get_cond_ids(self):
        return [entry['col_id'] for entry in self.dbclient["col_info"].find({}, {"_id": 0, "col_id": 1})]

    def get_corems(self):
        return [corem for corem in self.dbclient["corem"].find({}, {'_id': 1, "rows": 1, "corem_id": 1})]

    def no_col_resamples(self, col, nrows, nresamples):
        return self.dbclient.col_resample.find_one({"n_rows": nrows, "col_id": col, "resamples": {"$gte": nresamples}}) is None

    def find_gene_expressions(self, rows, cols):
        return pd.DataFrame(list(self.dbclient.gene_expression.find({"col_id": {"$in": cols}, "row_id": {"$in": rows}},
                                                                    { "_id": 0, "col_id": 1, "raw_expression": 1,
                                                                      "standardized_expression": 1})))

    def find_col_resamples(self, nrows, cols):
        return pd.DataFrame(list(self.dbclient.col_resample.find({"n_rows": nrows, "col_id": {"$in": cols}}, {"_id": 0})))

    def update_corem(self, corem, new_cols):
        self.dbclient['corem'].update({"_id": corem['_id']}, {"$set": {"cols": new_cols}})

    def get_row_maps(self):
        row2id = {}
        id2row = {}
        for i in self.dbclient.row_info.find({}, {"egrin2_row_name": "1", "row_id": "1"}):
            row2id[i["egrin2_row_name"]] = i["row_id"]
            id2row[i["row_id"]] = i["egrin2_row_name"]
        return row2id, id2row

    def get_column_maps(self):
        col2id = {}
        id2col = {}
        for i in self.dbclient.col_info.find({}, {"egrin2_col_name": "1", "col_id": "1"}):
            col2id[i["egrin2_col_name"]] = i["col_id"]
            id2col[i["col_id"]] = i["egrin2_col_name"]
        return col2id, id2col

    def num_row_co_occurence(self, rowname, row2id, id2row):
        """Given a row (gene), count all of the other rows that occur with it in a bicluster"""
        data = []
        for i in self.dbclient.bicluster_info.find({"rows": {"$all": [row2id[rowname]]}}, {"rows": "1"}):
            for j in i["rows"]:
                try:
                    data.append(id2row[j])
                except:
                    continue
        data_counts = pd.Series(data).value_counts()
        return data_counts

    def drop_row_rows(self):
        self.dbclient.row_row.drop()

    def update_row_row(self, keyrow_pk, subrow_pk, data_counts_norm, backbone_pval):
        self.dbclient.row_row.update({"row_ids": [subrow_pk, keyrow_pk]},
                                     {"$set": {"weight": data_counts_norm, "backbone_pval": backbone_pval}})

    def insert_row_row(self, rowrow_docs):
        self.dbclient.row_row.insert(rowrow_docs)

    def insert_corem(self, corem_docs):
        self.dbclient.corem.insert(corem_docs)

    def corem_sizes(self):
        return list(set([len(i["rows"] ) for i in self.dbclient["corem"].find({}, {"rows": 1})]))

