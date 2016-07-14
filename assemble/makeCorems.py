"""Generate corems using cMonkey ensemble MongoDB. Add them to an existing MongoDB"""
import os
import glob

import subprocess
import shutil
import decimal
from math import ceil
import logging
import numpy as np
import pandas as pd

from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def command_exists(command):
    return subprocess.call(["which", command], stdout=open(os.devnull, 'wb')) == 0

def check_c_code_exists():
    """checks whether the C commands required exist"""
    result = False
    if not command_exists("adjmat2wpairs"):
        logging.warn("You need to compile adjmat2wpairs.cpp to adjmat2wpairs and add its location to your path to detect corems")
        result = True

    if not command_exists("compute_tanimoto"):
        logging.warn("You need to compile compute_tanimoto.cpp to compute_tanimoto and add its location to your path to detect corems")
        result = True

    if not command_exists("cluster_communities"):
        logging.warn("You need to compile cluster_communities.cpp to cluster_communities and add its location to your path to detect corems")
        result = True

    if not command_exists("getting_communities"):
        logging.warn("You need to compile getting_communities.cpp to getting_communities and add its location to your path to detect corems")
        result = True

    return result


class MongoDB:
    """database interface to mongo"""
    def __init__(self, dbclient):
        self.dbclient = dbclient

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


class CoremMaker:

    def __init__(self, organism, db_client, backbone_pval, out_dir, n_subs, link_comm_score,
                 link_comm_increment, link_comm_density_score, corem_size_threshold,
                 n_resamples):

        self.organism = organism
        self.db_client = db_client
        self.row2id, self.id2row = db_client.get_row_maps()
        self.col2id, self.id2col = db_client.get_column_maps()
        self.backbone_pval = backbone_pval
        self.cFail = check_c_code_exists()

        self.out_dir = os.path.abspath(os.path.join(out_dir, "corem_data"))
        if not os.path.isdir(self.out_dir):
            os.makedirs(self.out_dir)
        logging.info("Corem data will be output to '%s'", self.out_dir)

        self.n_subs = n_subs
        self.link_comm_score = link_comm_score
        self.link_comm_increment = link_comm_increment
        self.link_comm_density_score = link_comm_density_score
        self.corem_size_threshold = corem_size_threshold
        self.cutoff = None
        self.n_resamples = n_resamples

    def __extract_backbone(self, data_counts):
        """Extract the significant elements from rBr co-occurrence matrix"""
        backbone_data_counts = data_counts.copy()

        def integrand(x, k):
             return np.power(1.0 - x, k - 2.0)

        for i in data_counts.index:
            k = len(data_counts)
            pval = 1 - (k - 1) * quad(integrand, 0, data_counts[i], args=(k))[0]
            backbone_data_counts[i] = pval
        return backbone_data_counts

    def __row_row(self):
        """Construct row-row co-occurrence matrix (ie gene-gene co-occurrence)"""

        # remove existing edgeList file if it exists
        if os.path.exists(os.path.abspath(os.path.join(self.out_dir, "edgeList"))):
            logging.info("Found edgeList file at '%s'. Removing it.",
                         os.path.abspath(os.path.join(self.out_dir, "edgeList")))
            os.remove(os.path.abspath(os.path.join(self.out_dir, "edgeList")))
            logging.info("Dropping row_row collection as a precaution.")
            self.db_client.drop_row_rows()

        def addToD(d, ind1, ind2, val):
            if ind1 not in d.keys():
                d[ind1] = {}
            d[ind1][ind2] = val
            return d

        def writeRowRow(f, row):
            """Only write rows with significant backbone pvals to edgeList file"""
            if float(row["backbone_pval"]) <= self.backbone_pval:
                f.write((" ").join([self.id2row[row["row_ids"][0]], self.id2row[row["row_ids"][1]], str(row["weight"]), "\n"]))

        def structureRowRow(key_row, sub_row, data_counts, data_counts_norm, backbone_pval):
            keyrow_pk = self.row2id[key_row]
            subrow_pk = self.row2id[sub_row]
            try:
                # check to see if this pair already exists and if current weight is greater
                weight = self.rowrow_ref[subrow_pk][keyrow_pk]

                if (data_counts_norm > weight) and (backbone_pval <= self.backbone_pval):
                    # if current weight is greater and backbone_pval is significant, update the weight in MongoDB
                    self.db_client.update_row_row(keyrow_pk, subrow_pk, data_counts_norm, backbone_pval)

                    self.rowrow_ref = addToD(self.rowrow_ref, subrow_pk, keyrow_pk, data_counts_norm)
                d = None
            except Exception:
                self.rowrow_ref = addToD(self.rowrow_ref, subrow_pk, keyrow_pk, data_counts_norm)
                d = {
                    "row_ids": [keyrow_pk, subrow_pk],
                    "counts": data_counts,
                    "weight": data_counts_norm,
                    "backbone_pval": backbone_pval
                }
            return d

        counter = 1

        # make a dictionary to keep track of rows in db and their weights
        # was too slow using mongoDB lookups...
        self.rowrow_ref = {}
        logging.info("Constructing row-row co-occurrence matrix. This will take some time...")

        for i in self.row2id.keys():
            if counter % 250 == 0:
                logging.info("%.2f percent done", round(float(counter) / len(self.row2id.keys()), 2) * 100.0)

            # check if already exists in DB
            data_counts = self.db_client.num_row_co_occurence(i, self.row2id, self.id2row)

            # set self counts to 0 and normalize other counts
            data_counts[i] = 0
            data_counts = data_counts[data_counts > 0]
            data_counts_norm = data_counts.copy()
            data_counts_norm = data_counts_norm / sum(data_counts_norm)

            # only keep values > 0
            backbone_data_counts = self.__extract_backbone(data_counts_norm)

            to_write = [structureRowRow(i, j, data_counts[j], data_counts_norm[j], backbone_data_counts[j])
                        for j in data_counts.index]
            to_write = [i for i in to_write if i is not None]

            # write edgeList file
            with open(os.path.abspath(os.path.join(self.out_dir, "edgeList")), mode="a+") as f:
                for j in to_write:
                    writeRowRow(f, j)

            self.db_client.insert_row_row(to_write)
            counter += 1

        # clean up
        del self.rowrow_ref

    def __run_corem_cscripts(self):

        def drange(start, stop = None, step = 1, precision = None):
            """drange generates a set of Decimal values over the
            range [start, stop) with step size step

            drange([start,] stop, [step [,precision]])"""

            try:
                _xrange = xrange
            except NameError:
                _xrange = range

            if stop is None:
                for x in _xrange(int(ceil(start))):
                    yield x
            else:
                # find precision
                if precision is not None:
                    decimal.getcontext().prec = precision
                # convert values to decimals
                start = decimal.Decimal(start)
                stop = decimal.Decimal(stop)
                step = decimal.Decimal(step)

                # create a generator expression for the index values
                indices = (i for i in _xrange(0, int(((stop-start) / step).to_integral_value())))

                # yield results
                for i in indices:
                    yield float(start + step*i)

        if self.cFail:
            logging.info("""Cannot detect corems because one or more community detection C++ scrips are either
(1) not compiled, (2) not in the $PATH variable, or (3) incorrectly named. Resolve previous warning.""")
            return None
        else:
            p = subprocess.Popen(["adjmat2wpairs", "edgeList", "0", "0"], cwd=os.path.abspath(self.out_dir))
            p.wait()

            ranges = list(range(0, len(self.row2id) + 1, int((len(self.row2id) + 1) / self.n_subs)))

            # make sure last is # genes
            ranges[len(ranges) - 1] = len(self.row2id) + 1

            # "0" parameter is for link similarity of Ahn et al. Use "1" for measure proposed by Kalinka et al
            command_template = ['compute_tanimoto', 'edgeList', str(self.link_comm_score), '0', len(self.row2id)]
            commands = []

            for i in range(0, len(ranges) - 1):
                cmd = list(command_template)
                cmd[3] = str(ranges[i])
                cmd[4] = str(ranges[i + 1])
                commands.append(cmd)

            # run in parallel
            processes = [subprocess.Popen(cmd, cwd=os.path.abspath(self.out_dir)) for cmd in commands]

            # wait for completion
            for p in processes:
                p.wait()

            # clean up
            with open(os.path.join(os.path.abspath(self.out_dir), "edgeList.tanimoto"), 'w' ) as outfile:
              for infile in glob.glob(os.path.join(os.path.abspath(self.out_dir), "edgeList.tanimoto_*")):
                shutil.copyfileobj(open(infile), outfile)

            p = subprocess.Popen(["rm edgeList.tanimoto_*"], cwd=os.path.abspath(self.out_dir), shell=True)
            p.wait()

            logging.info("Clustering link communities across thresholds defined by increment: %s", str(self.link_comm_increment))
            command_template = ['cluster_communities', 'edgeList', 1]
            commands = []

            for i in list(drange(self.link_comm_increment, 1 + self.link_comm_increment, self.link_comm_increment, 1)):
                cmd = list(command_template)
                cmd[2] = str(i)
                commands.append(cmd)

            # run in parallel
            count_1 = 0
            while count_1 < len(commands):
                # only use define number of subprocesses at a time
                count_2 = count_1 + self.n_subs
                processes = [subprocess.Popen(cmd, cwd=os.path.abspath(self.out_dir)) for cmd in commands[count_1: count_2]]

                # wait for completion
                for p in processes:
                    p.wait()

                count_1 = count_2
                if count_2 < len(commands):
                    logging.info("%.2f percent done", round(float(count_2) / len(commands), 2) * 100)

            # clean up
            with open(os.path.join(os.path.abspath(self.out_dir), "edgeList.density"), 'w') as outfile:
              for infile in glob.glob(os.path.join(os.path.abspath(self.out_dir), "edgeList.density_*")):
                shutil.copyfileobj(open(infile), outfile)

            p = subprocess.Popen(["rm edgeList.density_*"], cwd=os.path.abspath(self.out_dir), shell=True)
            p.wait()

            # choose cutoff
            density = pd.read_csv(os.path.join(os.path.abspath(self.out_dir), "edgeList.density"), sep="\t", header=None)

            # map density score to column
            score_map = {1: 2, 2: 4, 5: 6, 3: 7, 4: 8}

            vals = density[score_map[self.link_comm_density_score]]
            maxT = max(vals)
            maxT_ind = round(density.iloc[:, 0][vals.idxmax()], 4)

            try:
                # plot density scores
                black = [2, 4, 7]
                blue = [8, 6]
                lty = ['-', '--', '-.', ':']

                fig, ax1 = plt.subplots()
                l1 = ax1.plot(density.loc[:, 0], density.loc[:, black[0]], color="black", ls = lty[0], marker=".", ms=10, lw=2, label="(1) Unweighted (all)")
                l2 = ax1.plot(density.loc[:, 0], density.loc[:, black[1]], color="black", ls = lty[1], marker=".", ms=10, lw=2, label="(2) Unweighted (n>2)")
                l3 = ax1.plot(density.loc[:, 0], density.loc[:, black[2]], color="black", ls = lty[2], marker=".", ms=10, lw=2, label="(3) Weighted (local)")

                if score_map[self.link_comm_density_score] in black:
                    yl = ax1.plot(density.loc[:, 0], density.loc[:, score_map[self.link_comm_density_score]], color="red", ms=10, ls='-', marker=".", lw=2, label="Your choice")
                    ax1.axvline(maxT_ind, color='red', linestyle=':')

                ax1.set_xlabel('Threshold')
                # Make the y-axis label and tick labels match the line color.
                ax1.set_ylabel('Similarity Score', color='black')
                ax1.set_xlim((0, 1))
                for tl in ax1.get_yticklabels():
                  tl.set_color('black')

                ax2 = ax1.twinx()
                l4 = ax2.plot(density.loc[:, 0], density.loc[:, blue[0]], color="blue", ls=lty[0], marker=".", ms=10, lw=2, label= "(4) Weighted (global)")
                l5 = ax2.plot(density.loc[:, 0], density.loc[:, blue[1]], color="blue", ls=lty[1], marker=".", ms=10, lw=2, label="(5) Weighted (both)")

                if score_map[ self.link_comm_density_score ] in blue:
                    yl = ax2.plot(density.loc[:, 0 ], density.loc[:, score_map[self.link_comm_density_score]],
                                  color="red", ls='-', marker=".", ms=10, lw=2, label="Your choice")
                    ax2.axvline(maxT_ind, color='red', linestyle=':')

                ax2.set_ylabel('Similarity Score', color='blue')

                for tl in ax2.get_yticklabels():
                    tl.set_color('blue')

                #plt.plot( (maxT, maxT), (-100, 100), 'r--' )

                plt.legend(handles = [l1[0], l2[0], l3[0], l4[0], l5[0], yl[0]], loc = 2)
                plt.title('Corem community density at different similarity cutoffs')
                plt.tight_layout()
                pp = PdfPages(os.path.join(os.path.abspath(self.out_dir),"density_stats.pdf"))

                plt.savefig(pp, format='pdf')
                pp.close()

                logging.info("Threshold density plots written to: %s",
                             os.path.join(os.path.abspath(self.out_dir), "density_stats.pdf"))

            except Exception:
                logging.error("Could not produce density plot. Check matplotlib.")

            logging.info("Threshold is %.5f at cutoff = %.5f", round(maxT, 5), maxT_ind)
            self.cutoff = maxT_ind

            # end plot

    def __get_corems(self, cutoff=None):
        """load clusters at selected density"""
        if self.cutoff is None:
            if cutoff is None:
                return "Please provide a cutoff value"
            else:
                self.cutoff = cutoff

        if self.cFail:
            logging.error("""Cannot detect corems because one or more community detection C++ scrips are either
(1) not compiled, (2) not in the $PATH variable, or (3) incorrectly named. Resolve previous warning.""")
            return None

        # get communities at selected cutoff
        p = subprocess.Popen(["getting_communities", "edgeList", str(self.cutoff)], cwd=os.path.abspath(self.out_dir))
        p.wait()

        clusters = pd.read_csv(os.path.join(os.path.abspath(self.out_dir), "edgeList.communities_" + str(self.cutoff)),
                               sep="\t", header=None)
        clusters.columns = ["Gene1", "Gene2", "Community_ID", "Community_Density", "Community_Weighted_Density"]

        # find and keep communities with > 3 sig edges
        edgeCounts = clusters.Community_ID.value_counts()
        sig = edgeCounts[edgeCounts >= 3].index.tolist()
        sigClusters = clusters[clusters["Community_ID"].isin(sig)]

        # keep communities with > 0 density
        sigClusters = sigClusters[sigClusters["Community_Density"] > 0]

        # sort by density
        sigClusters = sigClusters.sort_values(['Community_Density','Community_Weighted_Density','Community_ID'],ascending=False)

        # rename corems
        clusterNameD = dict(zip(sigClusters.Community_ID.unique(), range(1, len(sigClusters.Community_ID.unique()) + 1)))
        sigClusters.Community_ID = [clusterNameD[i] for i in sigClusters.Community_ID]

        sigClusters.to_csv(os.path.join(os.path.abspath(self.out_dir), "edgeList.communities_" + str(self.cutoff) + "_FINAL.txt"),
                           sep="\t", index=False)
        return None

    def __add_corems(self):
        """Add corems to MongoDB. Will Only run if self.cutoff has been set by running C++ codes"""
        pd.options.mode.chained_assignment = None

        def coremStruct(corem, table):
            """MongoDB corem template"""
            logging.info("%d of %d corems completed", corem, len(table.Community_ID.unique()))

            sub_m = table.loc[table.Community_ID==corem, :]
            # translate names
            sub_m.loc[:, "Gene1"] = [str(self.row2id[i]) for i in sub_m.Gene1]
            sub_m.loc[:, "Gene2"] = [str(self.row2id[i]) for i in sub_m.Gene2]

            rows = list(set(sub_m.Gene1.unique().tolist() + sub_m.Gene2.unique().tolist()))
            rows = [int(i) for i in rows]
            rows.sort()

            edges = list(sub_m.Gene1  + "-" + sub_m.Gene2)

            d = {
                "corem_id": corem,
                "rows": rows,
                "cols": [],
                "edges": edges,
                "density": sub_m.Community_Density.unique()[0],
                "weighted_density": sub_m.Community_Weighted_Density.unique()[0]
            }

            return d

        if self.cutoff is None:
            return "Please run C++ code to determine a cutoff."

        corems = pd.read_csv(os.path.join(os.path.abspath(self.out_dir),
                                          "edgeList.communities_" + str(self.cutoff) + "_FINAL.txt"),
                             sep="\t")
        logging.info("Adding basic corem information to database")

        to_write = [coremStruct(i, corems) for i in corems.Community_ID.unique()]
        self.db_client.insert_corem(to_write)

    def make_corems(self):
        """top-level corem making function"""
        self.__row_row()
        self.__run_corem_cscripts()
        self.__get_corems()
        self.__add_corems()
