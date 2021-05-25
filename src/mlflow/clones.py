from src.utils.protocols_class import Protocol, Analysis, Var, helper_add_prefix_to_recursive_dict
import src.mlflow.mlflow_utils as mu

import seaborn as sns
import matplotlib.pyplot as plt
import copy
from mplh import cluster_help as ch
from mplh.fig_utils import helper_save as hs
from os.path import join
import pandas as pd
from src.external.pyvenn import venn
#from src.calculate_AF_by_cell import calculate_af
import numpy as np
import src.utils.variant_utils as vu
import click
import mlflow
import tempfile

import pickle

####
## dfs i/o
####
def read_af(files, sep='\t', to_transpose=False, **pdread_kwargs):
    """ Reads in and concat dfs from files dictionary

    :param files:
    :param sep:
    :return:
    """
    dfs = {}
    if type(files) != dict and type(files) != pd.Series:
        files = {f: f for f in files}
    for f in files:
        print(f)
        dfs[f] = pd.read_csv(files[f], sep=sep, **pdread_kwargs)
        if to_transpose:
            dfs[f] = dfs[f].transpose()
    return dfs


def combine_dfs(dfs, use_key=True, key_label="Key", naval=0,
                update_col_names=False, col_inds=None, axis=1):
    """

    :param dfs: dict of dfs, where keys are names to be used for a new column,
                and/or for add that as a suffix to the column names.
    """
    for key in dfs:
        if use_key:
            dfs[key][key_label] = key
        if update_col_names:
            dfs[key].columns = dfs[key].columns + "_" + key
        if col_inds is not None:
            dfs[key] = dfs[key].loc[:, dfs[key].columns.isin(col_inds)]
    return pd.concat(dfs.values(), axis=axis, sort=False).fillna(
        naval)  # return  pd.concat(dfs.values()).fillna(fillna)


####
## dfs transform+filtering
####
def drop_zero(df, axis=0):
    if axis == 0:
        return df.loc[(df == 0).all(axis=0)]
    else:
        return df[(df == 0).all(axis=1)]


def binarize(df, het_thresh):
    return df>het_thresh


def filter_variants(df, variant_inds):
    return df.loc[variant_inds]


def filter_variants_sums(df, variant_thresh, proportion_to_pass,
                         rm_high=True):
    high_variant_bool = ((df > variant_thresh).sum(axis=1) / df.shape[
        1] > proportion_to_pass)
    if rm_high:
        print('removing variants of high heteroplasmy')
        #print(f"Removing variants with greater than {variant_thresh} in {proportion_to_pass} frac of group")
        df = df.loc[~(high_variant_bool)]
    else:
        #print(
        #    f"Removing variants with less than {variant_thresh} in {proportion_to_pass} frac of group")
        df = df.loc[(high_variant_bool)]

    print(f"Number of variants above: {high_variant_bool.sum()}")
    print(f"Number of variants below: {(~high_variant_bool).sum()}")
    return df


def filter_cells(df, cell_inds):
    return


####
## Variants
####
def plot_variant_types(variant_df, out_f=None):
    f, ax = plt.subplots(2, 1)
    sns.countplot(variant_df["variant type"], ax=ax[0])
    sns.countplot(variant_df["variant change"], ax=ax[1])
    hs(out_f)
    return


def plot_variant_types_cat(variant_df, afs , out_f=None):
    sns.catplot(x="variant type", y="Total coverage",
                data=pd.concat((variant_df, afs.sum(axis=1)),sort=False,
                               axis=1).rename({0: "Total coverage"},
                                              axis=1))
    if out_f is not None:
        hs(out_f + ".variantType")
    sns.catplot(x="variant change", y="Total coverage",
                data=pd.concat((variant_df, afs.sum(axis=1)),sort=False,
                               axis=1).rename({0: "Total coverage"},
                                          axis=1))
    if out_f is not None:
        hs(out_f+".variantChange")


def get_variants_df(df, to_plot=False, out_f=None):
    variant_df = vu.type_of_variants(df.index.values)
    if to_plot:
        plot_variant_types(variant_df, out_f)
    return variant_df


def plot_overlap_variants(df=None, samples=None, variant_dict=None, out_f=None,
                          key="Condition",
                          variant_thresh=0.0, proportion_to_pass=1e-10):
    """ Plots venn diagram of overlapping samples.

    a. Does this by grouping the df by samples inds, and then getting which variants
    have all 0's in the condition.
    b. variant_dict where the keys are the samples and values are their variants
    """
    if variant_dict is not None:
        print("using variant dict")
    elif samples is not None and df is not None:
        variant_dict = {}
        for ind, val in samples.groupby(key):
            curr_df = df.loc[:, df.columns.isin(val.index)]
            variant_dict[ind] = filter_variants_sums(curr_df,
                                                     variant_thresh=variant_thresh,
                                                     proportion_to_pass=proportion_to_pass,
                                                     rm_high=False).index
    else:
        print("No list. Not plotting")
        return
    labels = venn.get_labels(list(variant_dict.values()))
    f, ax = venn.venn(labels, names=list(variant_dict.keys()), )
    # venn2([set(bin_J2_AF_by_cell.index.values),set(bin_P2_AF_by_cell.index.values)],set_labels=["P2","J2"])
    plt.title("Overlap of called variants")
    #f.show()
    hs(out_f)
    return


def plot_variants(af):
    return


####
## Plotting
####

def plot_depth_clustermap():
    return


def plot_variant_violinplots(df, variant_inds=None, cell_subsets=None,
                             use_joy=False):
    return


def plot_af_clustermap(df, vmax=1, f_save=None):
    g = ch.plot_cluster(df=df, row_meta=None, col_meta=None, fsave=None,
                        to_z=False, to_col_clust=True,
                        to_row_clust=True, name=None, col_names=True,
                        vmin=0, vmax=vmax, row_names=True,
                        to_legend=True, method="average",
                        white_name=None, metric='jaccard')
    plt.title("PBMC P cells-by-variants")
    hs(f_save)
    return g


def plot_heteroplasmy_shifts(afs, samples, samples_key="Condition",
                             out_f=None):
    """ Plot the heteroplasmy shift across samples

    :param afs:
    :param samples:
    :return:
    """

    var_het = pd.DataFrame(columns=["mu AF", "Condition"])
    for ind, val in samples.groupby(samples_key):
        curr_afs = afs.loc[:, val.index]
        var_mean = curr_afs.mean(axis=1)
        #pd.concat(var_het)
        var_het = pd.concat((var_het, pd.DataFrame({"mu AF": var_mean,
                                                    "Condition": ind,
                                                    "Var": var_mean.index})),
                            sort=False, axis=0)
        #print('var_het')
        #print(var_het)
        #var_het.loc[ind, "mu AF"] = var_mean
        #var_d[ind] = var_mean

    f, ax = plt.subplots()
    count = 0
    for ind, val in var_het.groupby("Var"):
        sns.lineplot(x="Condition", y="mu AF", data=val, ax=ax, ci=None,
                     color='blue')
        count += 1

    # sns.regplot(x="Condition", y="mu AF", data=var_het, x_jitter=.1)
    hs(out_f)
    return var_het

####
## cell indices subgroups for Donors and clones
####
def separate_donors(df, donor_df, remove_variants=True,
                    variant_thresh=0.4, donor_thresh=0.9):

    donors = {}
    metrics = {}
    for d, val in donor_df.groupby("Donor"):
        print(f"Donor {d}")
        curr_df = df.loc[:, val.index]
        if remove_variants:
            print(
                f"Removing variants with greater than {variant_thresh} in {donor_thresh} frac of group")
            high_variant_bool = (
                        (curr_df > variant_thresh).sum(axis=1) /
                        curr_df.shape[1] > donor_thresh)
            high_vars = curr_df.loc[high_variant_bool].index
            print(f"Number of variants removed: {len(high_vars)}")
            curr_df = curr_df.loc[~(high_variant_bool)]
        donors[d] = curr_df
        metrics[f"high_variants_removed_donor{d}"] = len(high_vars)
    return donors, metrics


def load_lineage(in_files_dict, lin_col="Lineage"):
    """Load the cell indices for each cluster and each sample

    :param in_files_dict: Keys are names of cluster, values are the file paths
    :return dictionary where keys are names of cluster and values is list
    of cell names.

    THIS taken from Analysis/lineage_and_peakclusters/compare_lineage_clusters.py
    """
    lineages = {}
    for i in in_files_dict:
        lineages[i] = pd.read_csv(in_files_dict[i])[["ID"]]
        lineages[i][lin_col] = i
    return pd.concat(lineages, ignore_index=True)


def get_lineage_files(in_dir, is_donor=True, N_DONOR=None):
    return []


def d_from_lst(keys, vals):
    d ={val:vals[ind] for ind, val in enumerate(keys)}
    return d


def prep_files(f_lists, samples):
    return [d_from_lst(samples, x) for x in f_lists]


def load_donors(donor_dir, N_DONORS, out_f=None):
    pairs = dict()
    for d in range(N_DONORS):
        # add to the pair dictionary the {donor: file}
        pairs[d] = join(donor_dir,
                         f"cell_labels.donor{d}.txt")

    donors_df = load_lineage(pairs, lin_col="Donor")
    donors_df["Condition"] = donors_df["ID"].apply(lambda x: x.split("_")[-1])
    donors_df = donors_df.set_index("ID")
    sns.countplot(data=donors_df, x="Donor", hue="Condition")
    hs(out_f)
    return donors_df



####
## Clone analysis
####
def call_clones(df, how='graphcluster'):
    """
    how: {'graphcluster', 'variant', 'variant_noDuplicates' }
    no_dup: {None, "smallest_group", "max_af", "smallest_af"}
    """
    return


def plot_mitomaps(af, cov, samples, out_f="",
                  bin_af=None, max_to_plot=1000):
    """ Plot the AF & heatmap coverages

    :param af:
    :param cov:
    :param samples:
    :param out_f:
    :param bin_af:
    :param max_to_plot:
    :return:
    """
    to_save = True
    if out_f == "":
        to_save = False
    if max_to_plot<af.shape[1]:
        inds = af.sample(n=max_to_plot, axis=1).columns
        af = af.loc[:, inds]
        cov = cov.loc[:, inds]
        if bin_af is not None:
            bin_af = bin_af.loc[:, inds]
    if max_to_plot<af.shape[0]:
        rows = af.sample(n=max_to_plot, axis=0).index
        af = af.loc[rows]
        cov = cov.loc[rows]
        if bin_af is not None:
            bin_af = bin_af.loc[:, rows]

    g = ch.plot_cluster(df=af, fsave=out_f + ".af.png", to_z=False,
                        to_col_clust=True, col_meta=samples,
                        to_row_clust=True, name="AF Pos-by-cell",
                        col_names=False, row_names=False,
                        to_legend=True, method="average",
                        white_name=None, metric='euclidean',
                        to_save=to_save)

    # Use same positions as the af
    ch.plot_cluster(df=cov.iloc[
        g.dendrogram_row.reordered_ind, g.dendrogram_col.reordered_ind],
                    fsave=out_f + ".cov.png", to_z=False,
                    to_col_clust=False, col_meta=samples,
                    to_row_clust=False, name="Depth Pos-by-cell",
                    col_names=False, row_names=False, to_legend=True,
                    to_save=to_save)

    if bin_af is not None:
        # Use same positions as the af
        ch.plot_cluster(figsize=(7, 4), dpi=300, df=bin_af.iloc[
            g.dendrogram_row.reordered_ind, g.dendrogram_col.reordered_ind],
                        fsave=out_f + ".afBinary.png", to_z=False,
                        to_col_clust=False, col_meta=samples,
                        to_row_clust=False, name="Depth Pos-by-cell",
                        col_names=False, row_names=False,
                        to_legend=True, method="average",
                        white_name=None, to_save=to_save)

    #sns.distplot(cov.sum(axis=0))
    #sns.distplot(np.log10(cov.sum(axifs=0)+1))
    # sns.displot(data=pd.merge(np.log2(cov.sum(axis=0) + 1), samples,
    #                           how='inner'), x=0, hue="Condition",
    #             kind="kde")
    return g




#######################################################################
## Variable classes
#######################################################################
class Cov(Var):
    def __init__(self, file=None, var=None):
        super().__init__(file, var)
        return

    def load(self, file=None, return_var=True, type='filters'):
        if file is not None:
            self.file = file
        if self.file is None:
            raise ValueError("Please set file first")

        if type=="af_by_cell":
            self.var = read_af(self.file)
        elif type=="filters":
            #The input is formatted slightly differently
            self.var = read_af(self.file, index_col=0, to_transpose=True)
        if return_var:
            return self.var
        return



#######################################################################
## Protocol classes
#######################################################################
class MTExp(Protocol):
    def __init__(self, config_f=None, config=None,
                 pipeline_id=-1,
                 input_d=None,
                 output_d=None):
        super().__init__(config_f, config,
                 pipeline_id,
                 input_d,
                 output_d)
        self.analysis = MTPipeline
        return

    def run(self):
        """ Method called in the entry points to run protocol.
        :return:
        """
        self.initialize()
        #self.load_input_files()
        self.init_mlflow(mlflow_tracking_uri="")
        if self.parameters['input_mode'] == "mgatk_vireoDonors":
            self.run_mgatk_vireoDonors()
            is_completed=True
            # try:
            #     self.run_mgatk_vireoDonors()
            #     is_completed=True
            # except:
            #     is_completed=False
        else:
            print("not a proper input mode", self.parameters['input_mode'])
            return
        if is_completed:
            self.is_completed=True
        else:
            print("Protocol failed")
        #self.save(vars=False)
        return


    def mlflow_log(self):
        """
        Logs all the mlflow artifacts. This is called by run, but can
        have different logs for different steps.
        :return:
        """
        return


    @staticmethod
    def input_keys():
        return {"in_af_d": "af",
                "in_mgatk_d": "mgatk_af",
                "in_cov_d": "cov",
                "donor_dir": "donors_df"}

    @staticmethod
    def var_loader(var, fname, params=None):
        if var == "af":
            return read_af(fname)
        elif var == "mgatk_af":
            return read_af(fname)
        elif var == "cov":
            print('cov fname', fname)
            return Cov(file=fname).load() # mttrace.load_output_files(fname, var_keys=["cov"])
        elif var == "donors_df":
            return Donors.load_output_files({"donor_dir": fname},
                                            var_keys=["donor_dir"],
                                            params={"N_DONORS": params["N_DONORS"]})
        else:
            return f"No loader for variable {var}"
        # cov = read_af(in_cov_d)


    @staticmethod
    def output_files(outdir="", files=None, params=None, name=""):
        #ndonor = params["N_DONORS"]
        return {"results": {"overlapvar":"overlap_variants.txt",
                            "d_overlapvar":"donors/overlap_variants",
                            "d_clone_labels": "d_clone_labels.txt",
                            "d_clone_variants": "d_clone_variants.txt"},
                "data": {},
                "plots": {"varmeta":"variantMeta.png",
                          "varmetacounts":"variantMetaCounts.png",
                          "mtheat":"MT_heatmap.png",
                          "overlap":"overlapVariants.png",
                          "d_mtheat": "donors/MT_heatmap",
                          "d_varmeta":"donors/variantMeta",
                          "d_varmetacounts": "donors/variantMetaCounts",
                          "d_overlap": "donors/overlapVariants",
                          "d_hetshift": "donors/hetshifts"}}


    def test_output_files(self):
        return


    #######
    ## ENTRY: Preprocess
    #######
    def preprocess(self, mode):
        # Load MTExp info
        self.save_start()
        if mode == "filters":
            self.mgatk_preprocess()
        else:
            raise ValueError(f"{mode} not implemented for preprocess")
    def mgatk_preprocess(self):
        print("loading mgatk variables")
        self.load_mgatk()
        af = self.vars['af']
        cov = self.vars['cov']
        mgatk_af = self.vars['mgatk_af']

        ## Allele frequency processing
        # Filter for mgatk_afs
        mgatk_afs = {}
        variant_d = {}
        for s in af:
            mgatk_afs[s] = filter_variants(af[s], mgatk_af[s].index)
            variant_d[s] = mgatk_af[s].index.values
        self.vars['variant_d'] = variant_d
        # Drop any cells with no counts, which shouldnt really happen
        for s in mgatk_af:
            mgatk_af[s] = drop_zero(mgatk_af[s], axis=1)

        # Merge the AFs and add the sample name to the cell ID suffix
        self.vars["merged_afs"] = combine_dfs(copy.deepcopy(mgatk_afs),
                                 use_key=False, fillna=0,
                                 update_col_names=True)

        # Merge coverage and extract only the variant positions
        merged_covs = combine_dfs(copy.deepcopy(cov), use_key=False,
                                  fillna=0, update_col_names=True,
                                  col_inds=self.vars["merged_afs"].columns)
        self.vars["bin_af"] = binarize(self.vars["merged_afs"],
                                       het_thresh=0.001)
        self.vars["merged_covs"] = np.log2(merged_covs + 1)
        self.save_intermediate()
        return


    def load_mgatk(self):
        self.load_input_files(
            var_keys=["in_af_d", "in_mgatk_d", "in_cov_d"], return_vars=False)
        return

    @staticmethod
    def intermediate_files(files_dict=None):
        if files_dict is None:
            return {"merged_afs": "merged_afs.csv",
                    "merged_covs": "merged_covs.csv",
                    "bin_af": "bin_af.csv",
                    "variant_d": "variant_d.p"
                    }
        return {}

    def save_start(self):
        """ Save mlflow to start
        :return:
        """
        #with mlflow.start_run(nested=True):
        #out_dir = tempfile.mkdtemp()
        with tempfile.TemporaryDirectory() as td:
            with open(join(td, "start.txt"),'w') as f:
                f.write(self.config_f)
            mlflow.log_artifacts(local_dir=td,
                                 artifact_path="data/start")

    def save_intermediate(self):
        files = self.intermediate_files()
        #with mlflow.start_run(nested=True):
            #out_dir = tempfile.mkdtemp()
        with tempfile.TemporaryDirectory() as td:
            for var in files:
                if files[var].split(".")[-1] == "p":
                    pickle.dump(obj=self.vars[var],
                                file=open(join(td,files[var]),'wb'))
                else:
                    self.vars[var].to_csv(join(td,files[var]))
            mlflow.log_artifacts(local_dir=td,
                                 artifact_path="data/intermediate")
        return

    def load_intermediate(self, d="."):
        files = self.intermediate_files()
        print('files', files)
        print('d', d)
        for var in files:
            if files[var].split(".")[-1] == "p":
                self.vars[var] = pickle.load(file=open(join(d,files[var]), 'rb'))
            else:
                self.vars[var] = pd.read_csv(join(d,files[var]), index_col=0)
        return


    #######
    ## ENTRY 2: Results
    #######
    def results(self, data_uri):
        # Load donors
        self.load_input_files(var_keys=["donor_dir"],
                              return_vars=False)
        donors_df = self.vars['donors_df']['donors_df']


        self.load_intermediate(d=data_uri)
        merged_afs = self.vars["merged_afs"]
        merged_covs = self.vars["merged_covs"]
        #bin_af = self.vars["bin_af"]
        variant_d = self.vars['variant_d']


        # Create variant meta df with Position, Ref, Alt
        # And plot
        variant_df = get_variants_df(merged_afs, to_plot=False)
        # plot the heatmaps and variant info
        self.plot_mtexp(merged_afs, merged_covs, donors_df,
            variant_df, variant_d, "figures")

        ####
        ## Multiplex Donor extraction
        ####
        # Default parameters here
        donors, don_metrics = separate_donors(merged_afs, donors_df,
                                 remove_variants=True,
                                 variant_thresh=0.4, donor_thresh=0.75)
        thresholds = [(0.4,0.5), (0.1, 0.9), (0.4,0.9)]
        for t in thresholds:
            with mlflow.start_run(nested=True) as child_run:
                _ , metrics = separate_donors(merged_afs, donors_df,
                                    remove_variants=True,
                                    variant_thresh=t[0],
                                    donor_thresh=t[1])
                mlflow.log_metrics(metrics)

            ####
            ## Plotting Donors
            ####
            # Plot MT maps, including AF,
            # Plot variant overlap and heteroplasmy shifts between conditions
            for d in donors:
                #with tempfile.TemporaryDirectory() as out_dir:
                curr_af = donors[d]
                curr_cov = merged_covs.loc[curr_af.index, curr_af.columns]
                curr_donors_df = donors_df.loc[curr_af.columns]
                curr_var_df = variant_df.loc[curr_af.index]
                self.plot_mtexp(curr_af, curr_cov, curr_donors_df,
                                curr_var_df, None, f"figures/donors/donor{d}")
        return

    def plot_mtexp(self, merged_afs, merged_covs, donors_df, variant_d,
                   variant_df, artifact_path):
        print('plot_mtxp')
        print('af',merged_afs.head())
        print('cov',merged_covs.head())
        print('donors_df', donors_df.head())
        with tempfile.TemporaryDirectory() as td:
            ###
            # Plotting MT Maps of merged conditions
            ###

            plot_mitomaps(merged_afs, merged_covs, samples=donors_df,
                          out_f=join(td, "MT_heatmap"))
            plot_variant_types_cat(variant_df, merged_covs,
                                   out_f=join(td, 'variantMetaCounts'))
            plot_overlap_variants(df=None, samples=None,
                                  variant_dict=variant_d,
                                  out_f=join(td, 'overlapVariants'),
                                  key="Condition")
            plot_heteroplasmy_shifts(merged_afs, donors_df,
                                     samples_key="Condition",
                                     out_f=join(td,
                                                f'HetShift'))

            mlflow.log_artifacts(td, artifact_path=artifact_path)
        return


class Donors(Protocol):
    def __init__(self):
        super().__init__()
        return

    @staticmethod
    def output_files(outdir="", files=None, params=None, name=""):
        out_d = {"donor_dir": outdir }
        return out_d

    @staticmethod
    def load_output_files(output_d, var_keys=None, params=None):
        vars = {}
        if var_keys is None:
            var_keys = list(output_d.keys())
        if params is None:
            params = {}
        for v in var_keys:
            if v == "donor_dir":
                print('output_d', output_d)
                vars["donors_df"] = load_donors(output_d[v],
                                                   N_DONORS=params["N_DONORS"])
        return vars



class mttrace(Protocol):
    def __init__(self):
        super().__init__()
        return

    @staticmethod
    def output_files(outdir="", files=None, params=None, name=""):
        if files is None:
            files = {}
        files["in_cov_d"] = {}
        files["in_mgatk_d"] = {}
        files["in_af_d"] = {}

        for f in files["cov_dir"]:
            files["in_cov_d"][f] = join(files["cov_dir"][f], "af_by_cell.DP.tsv")
            files["in_mgatk_d"][f] = join(files["cov_dir"][f],
                                        "filter_mgatk", f"{f}_lowC{params['lowC']}.af.mgatk.tsv")
            files["in_af_d"][f] = join(files["cov_dir"][f],
                                        "filter_mgatk", f"{f}_lowC{params['lowC']}.af.tsv")
        return helper_add_prefix_to_recursive_dict(d=files,
                                                   in_dir=outdir)


class mgatk(Protocol):
    def __init__(self):
        super().__init__()
        return
    @staticmethod
    def output_files(outdir="", files=None, params=None, name=""):
        if files is None:
            files = {}
        if params is not None and len(params) >0:
            p = params
            cov_d = f"{p['results']}/{p['sample']}/MT/cellr_{p['cellr_bc']}/{p['sample']}_{p['num_read']}/filters"
            cov_d = f"{cov_d}_minC{p['min_cells']}_minR{p['min_reads']}_topN{p['topN']}_hetT{p['het_thresh']}_hetC{p['min_het_cells']}_hetCount{p['het_count_thresh']}_bq{p['bq_thresh']}/filter_mgatk"
            if 'cov_d' not in files:
                files["cov_d"] = cov_d
        # return {"results": ["A.csv",
        #                     "B.csv"],
        #         "data": ["A.png", "B.png"],
        #         "plots": ["A1.png", "B2.png"]}
        if "in_mgatk_d" not in files:
            files["in_mgatk_d"] = outdir
        return files

    def test_output_files(self):
        return


#######################################################################
## Analysis class
#######################################################################
class MTPipeline(Analysis):
    def __init__(self,params_f, name, description=None):
        super().__init__(params_f, name, description)
        return


    @staticmethod
    def protocol_file_setup(outdir="", protocol=None, params=None,
                            files=None):
        "MGATK"
        if params is None:
            params = {}
        loaders = {"MTExp" : MTExp.output_files,
                   "mgatk" : mgatk.output_files,
                   "mttrace": mttrace.output_files,
                   "Donors": Donors.output_files
                   }

        if protocol is None:
            return loaders
        return loaders[protocol](outdir=outdir, params=params, files=files)
    @staticmethod
    def protocol_runner(protocol=None):
        runners = {"MTExp": MTExp().run()}
        if protocol is None:
            return runners
        return runners[protocol]


