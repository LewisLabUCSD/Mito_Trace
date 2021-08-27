import pickle
from vireoSNP import Vireo, vireo_wrap
from vireoSNP.plot.base_plot import heat_matrix
import vireoSNP
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mplh.cluster_help as ch
import os
from os.path import join, exists
from joblib import Parallel, delayed
from scipy.io import mmread
from src.utils.data_io import wrap_write_mtx_df
from icecream import ic


def run_single(ad_dp, n):
    _modelCA = Vireo(n_var=ad_dp[1].todense().shape[0],
                    n_cell=ad_dp[1].todense().shape[1], n_donor=n, n_GT=2,
                    fix_beta_sum=False, ASE_mode=True)
    _modelCA.set_prior(beta_mu_prior=np.array([[0.01, 0.5]]))
    _modelCA.fit(ad_dp[0], ad_dp[1], min_iter=20, verbose=False)
    return _modelCA


def run_vireo(ad, dp, k, out_f=None, plot_qc=False, n_cores=24,
              n_initials=50):
    _ELBO_mat = []
    _models_all = Parallel(n_jobs=n_cores)(delayed(run_single)((ad, dp), k) for i in range(n_initials))

    _elbo_temp = [x.ELBO_[-1] for x in _models_all]
    _idx = np.argmax([x.ELBO_[-1] for x in _models_all])
    modelCA = _models_all[_idx]
    _losses = modelCA.ELBO_

    print([x.ELBO_[-1] for x in _models_all])
    if plot_qc:
        plt.figure(figsize=(11, 4))
        plt.subplot(1, 2, 1)
        plt.hist([x.ELBO_[-1] for x in _models_all])
        plt.ylabel("Frequency")
        plt.xlabel("ELBO with multiple initializations")
        plt.subplot(1, 2, 2)
        plt.plot(_losses)
        plt.xlabel("Iterations")
        plt.ylabel("ELBO in a single run")
        plt.tight_layout()
        #plt.show()
    if out_f is not None:
        pickle.dump(modelCA, open(out_f + f".clones{k}.modelCA.p", "wb"))
        plt.savefig(out_f + ".png")
        #plt.close()
    return modelCA, _elbo_temp

#######################################################################
## 2 extractions- extract_clusters doesnt get the cell IDs,
## separate_donors does. Need to merge at some point, but for now
## separate_donors for the multiplex, extract_clusters for lineage
def extract_clusters(modelCA, prob_thresh, doublet_thresh, doublet_prob,
                     cells_meta,
                     outdir=None, out_f="lineage"):
    low_conf_cells = np.flatnonzero(doublet_prob > doublet_thresh)
    cell_clusters = dict()

    cells_meta[out_f] = np.nan
    cells_meta[f"{out_f} "] = np.nan
    for n in range(modelCA.ID_prob.shape[1]):
        # Drop low probability and/or high doublet probability
        cell_clusters[n] = np.flatnonzero(
            (modelCA.ID_prob[:, n] > prob_thresh))
        cell_clusters[n] = cell_clusters[n][
            ~(np.isin(cell_clusters[n], low_conf_cells))]

        if cells_meta is not None:
            curr_cells_meta = cells_meta.iloc[cell_clusters[n]].copy()
            curr_cells_meta = curr_cells_meta.reset_index()

        cell_clusters[n] += 1
        # Get their IDs
        if outdir != "" and exists(outdir):
            curr_out = join(outdir, f"{out_f}_donor{n}_cells.txt")
            curr_str = "\n".join(cell_clusters[n].astype(str))
            with open(curr_out, "w") as f:
                f.write(curr_str)
            if cells_meta is not None:
                curr_cells_meta.to_csv(
                    join(outdir, f"{out_f}_{n}.labels.txt"))
    return cell_clusters


def separate_donors(AD, DP, modelCA, cells_meta, outdir,
                    doublet_prob,
                    prob_thresh=0.9, doublet_thresh=0.9,
                    cells_ind_col='new index',
                    out_name="donor", prefix="",
                    cells_filt_col=None, cells_filt_val=None, vars_meta=None):
    """Separates the matrices and labels by donor using the multiplex output.

    :param AD: sparse matrix, position-by-cell.
    :param DP: sparse matrix, position-by-cell.
    :param modelCA:
    :param cells_meta:
    :param outdir:
    :param N_DONORS:
    :param doublet_prob:
    :param prob_thresh:
    :param doublet_thresh:
    :return:
    """
    low_conf_cells = np.flatnonzero(doublet_prob > doublet_thresh)
    cell_clusters = dict()
    ic("Before filtering")
    ic(cells_meta.shape)

    # filter for certain columns
    if cells_filt_col is not None and cells_filt_val is not None:
        cells_meta = cells_meta[cells_meta[cells_filt_col]==cells_filt_val]
    cells_meta = cells_meta.reset_index() # make 0-based
    ic("After filtering")
    ic(cells_meta.shape)
    allAD = pd.DataFrame(AD.todense())
    allDP = pd.DataFrame(DP.todense())
    ic(allAD.shape)
    #ic(allAD.head())

    #assert((min(cells_meta.index)==1))
    cells_meta[out_name] = np.nan
    cells_meta[f"{out_name}_index"] = np.nan # Add additional index
    # For each donor, extract their cells for their matrices and labels
    # print('cell meta')
    # print(cells_meta.head())
    for n in range(modelCA.ID_prob.shape[1]):
        # Drop low probability and high doublet probability
        cell_clusters[n] = np.sort(np.flatnonzero(
            (modelCA.ID_prob[:, n] > prob_thresh)))
        cell_clusters[n] = np.sort(cell_clusters[n][
            ~(np.isin(cell_clusters[n], low_conf_cells))])

        ## Extract specific AD+DP indices and re-index to 1-based
        # 1. make dense 0-based and take only the cluster indices
        curr_ad = allAD.loc[:, cell_clusters[n]].reset_index().melt(
            id_vars='index', var_name="Cell",
            value_name="Count").rename({"index": "Position"}, axis=1)

        curr_dp = allDP.loc[:, cell_clusters[n]].reset_index().melt(
            id_vars='index', var_name="Cell",
            value_name="Count").rename({"index": "Position"}, axis=1)

        #ic(curr_dp["Cell"].unique())

        # 2. Then drop indices with AD=0 to make sparse to remove positions/cells.
        #    This may remove some cells
        curr_dp = curr_dp.loc[~(curr_ad["Count"] == 0)]
        curr_ad = curr_ad.loc[~(curr_ad["Count"] == 0)]

        #i for ind, valcell_clusters[n]
        # 3. Assign cells_meta labels and indices based on curr_dp 0-based indices.

        # If any labels were removed from step 2, remove
        keep_inds = np.sort(list(set(cell_clusters[n]).intersection((set(curr_ad["Cell"].values)))))
        ic(len(curr_dp["Cell"].unique())-len(keep_inds))
        ic(cells_meta.shape)
        cells_meta.loc[keep_inds, out_name] = n # df is 1-based index but model is 0-based
        cells_meta.loc[keep_inds, f"{out_name}_index"] = np.arange(1,len(keep_inds)+1) # Set new index to 1-based
        cells_meta = cells_meta.astype({f"{out_name}_index": "Int64", out_name: "Int64"})


        # 4. Update the cell and position maps using dp, making it 1-based
        #curr_cell_map = {val-1: val for val in cells_meta[f"{out_name}_index"].values} #1-based to 0-based
        curr_cell_map = {val: ind + 1 for ind, val in
                         enumerate(np.sort(curr_dp["Cell"].unique()))}
        curr_pos_map = {val: ind + 1 for ind, val in enumerate(
            np.sort(curr_dp["Position"].unique()))}
        curr_ad["Cell"] = curr_ad["Cell"].map(curr_cell_map)
        curr_ad["Position"] = curr_ad["Position"].map(curr_pos_map)
        curr_dp["Cell"] = curr_dp["Cell"].map(curr_cell_map)
        curr_dp["Position"] = curr_dp["Position"].map(curr_pos_map)

        #ic(curr_cell_map)
        # 5. Get only the passed variants
        if vars_meta is not None:
            curr_vars_meta = vars_meta.iloc[list(curr_pos_map.keys())]

        print(f"{out_name} {n}: {len(cell_clusters[n])} cells ")
        print(curr_dp.shape)
        print(curr_ad.shape)

        # Get cell labels for specific lineages
        curr_cells_meta = cells_meta.loc[keep_inds].copy()
        curr_cells_meta = curr_cells_meta.reset_index(drop=True) #make it 0-based
        # Get their IDs
        if outdir != "" and exists(outdir):
            curr_out = join(outdir, f"{prefix}{out_name}{n}_cells.txt")
            curr_str = "\n".join(cell_clusters[n].astype(str))
            with open(curr_out, "w") as f:
                f.write(curr_str)
            curr_cells_meta.to_csv(
                join(outdir, f"{prefix}{out_name}{n}.labels.txt"))
            curr_cells_meta.to_csv(
                join(outdir, f"cell_labels.{prefix}{out_name}{n}.txt"))
            #print('curr_ad', curr_ad.head())
            wrap_write_mtx_df(outdir, curr_ad, curr_dp,
                              oth=None, to_rm=True, prefix=f"{prefix}{out_name}{n}",
                              columns=('Position', 'Cell', 'Count'))
            if vars_meta is not None:
                curr_vars_meta.to_csv(join(outdir, f"{prefix}{out_name}{n}.vcf"), sep='\t', index   =False)

    # cells_meta = cells_meta.sort_values(f"{out_name}_index")
    # cells_meta = cells_meta.set_index(
    #     cells_ind_col)  # Set to original col1-based index
    if outdir != "" and exists(outdir):
        cells_meta.to_csv(join(outdir, f"{prefix}cells_meta.tsv"),sep='\t', index=False)
        AF_SNPs = np.sum(
            modelCA.GT_prob * np.expand_dims(modelCA.beta_mu, 1),
            axis=2)
        pd.DataFrame(AF_SNPs, columns=[f"Cluster {x}" for x in
                                       range(AF_SNPs.shape[1])]).to_csv(
            join(outdir, "AF_SNPs.csv"), index=False)
    return cell_clusters


#######################################################################
def plot_vireo_out(modelCA, out_f, to_sqrt=False, labels=None,
                   doublet_prob=None):
    if labels is not None:
        clust_df = pd.DataFrame(
            modelCA.ID_prob + np.random.uniform(low=0.0, high=0.0000005,
                                                size=modelCA.ID_prob.shape),
            columns=[f"AF {x + 1}" for x in
                     np.arange(modelCA.ID_prob.shape[1])],
            index=labels.index)
        if doublet_prob is not None:
            clust_df = pd.concat((clust_df, pd.DataFrame(doublet_prob,
                                                         columns=[
                                                             "Doublet"])),
                                 axis=1)

        #curr = clust_df.sample(n=min(1000, clust_df.shape[0]))
        rand_df = clust_df.sample(n=min(1000, clust_df.shape[0]))

        rand_labels = labels.loc[rand_df.index]

        print(rand_labels.head())
        ch.plot_cluster(rand_df, cmap='Oranges', alpha=0.8,
                        to_row_clust=True, to_col_clust=False,
                        row_meta=rand_labels, to_legend=True,
                        white_name=None, yticklabels=False)
        plt.suptitle(f"Cell-cluster probability {clust_df.shape[0]}\nmax 1000 cells shown")
        plt.savefig(out_f + ".labels.png")  # plt.close()

    AF_SNPs = np.sum(
        modelCA.GT_prob * np.expand_dims(modelCA.beta_mu, 1), axis=2)
    print("AF_SNPs shape", AF_SNPs.shape)

    ch.plot_cluster(pd.DataFrame(AF_SNPs, columns=np.arange(1,AF_SNPs.shape[1]+1)), cmap='Blues', alpha=0.8,
                    to_row_clust=True, to_col_clust=True,
                    to_legend=True, white_name=None)
    plt.savefig(out_f + ".variants.labels.png")
    # plt.show()
    if doublet_prob is not None and labels is not None:
        return clust_df, AF_SNPs
    return None


def run_elbo(ad, dp, out_f, cells_meta,
             n_clone_list=None, save_clusters=True,
             labels=None, rerun_model=True, n_cores=16, n_initials=50):
    """ Sweep n-clones parameter in vireo and plot the ELBO

    :param ad: Dense allele depth matrix
    :param dp: Dense total depth matrix
    :param out_f: Output prefix to save, include directory
    :param cells_meta: Labels
    :param n_clone_list:
    :param save_clusters:
    :param labels:
    :param rerun_model:
    :param n_cores:
    :param n_initials:
    :return:
    """
    if n_clone_list is None:
        n_clone_list = np.arange(2, 5)

    _ELBO_mat = []

    for k in n_clone_list:
        print('lineages', k)
        curr_clone_model_f = out_f + f".clones{k}.modelCA.p"
        curr_lineage_f = out_f + f".clones{k}"
        if os.path.exists(curr_clone_model_f) and (not rerun_model):
            print('Already ran. Loading model')
            try:
                modelCA, _elbo_temp = pickle.load(open(curr_clone_model_f, "rb"))
            except TypeError:
                modelCA, _elbo_temp = run_vireo(ad, dp, k, out_f,
                                                n_cores=n_cores,
                                                n_initials=n_initials)
        else:
            modelCA, _elbo_temp = run_vireo(ad, dp, k, out_f,
                                n_cores=n_cores, n_initials=n_initials) #+ f".clones{k}.modelCA.p",
        # Add best loss to ELBO vector
        _ELBO_mat.append(_elbo_temp)
        print('_ELBO_mat', _ELBO_mat)
        # Run plot_vireo_out and extract the cell IDs as well
        if save_clusters and not (out_f == ""):
            print('saving lineage tree for file:', curr_lineage_f)
            try:
                doublet_prob = \
                modelCA.predict_doublet(ad, dp, update_GT=False,
                                        update_ID=False)[0].sum(axis=1)
            except AttributeError:  # New version of Vireo 2021
                doublet_prob = \
                vireoSNP.utils.vireo_doublet.predict_doublet(modelCA,
                                                             ad, dp,
                                                             update_GT=False,
                                                             update_ID=False)[
                    0].sum(axis=1)
            plot_vireo_out(modelCA, out_f=curr_lineage_f, to_sqrt=True,
                           labels=labels, doublet_prob=doublet_prob)
            separate_donors(ad, dp, modelCA, cells_meta, out_f,
                            doublet_prob, prob_thresh=0.9,
                            doublet_thresh=0.9,
                            cells_ind_col='donor index',
                            out_name="lineage")

    print('_ELBO_mat', _ELBO_mat)
    if len(_ELBO_mat) > 1:
        f = plt.figure()
        plt.boxplot(_ELBO_mat)
        plt.plot(np.arange(1, len(n_clone_list) + 1),
                 np.max(_ELBO_mat, axis=1))
        # plt.xticks(n_clone_list)
        plt.gca().set_xticklabels(n_clone_list)
        plt.ylabel("ELBO")
        plt.xlabel("n_clones")
        if out_f != "":
            plt.savefig(out_f+"_lineage_elbow.png")
        #plt.show()
    return _ELBO_mat, n_clone_list


def run_lineage(n, indir, outdir, n_clone_list=(5, 10, 20, 40, 100),
                rerun_model=False):
    """ Lineage trace on a donor across various clone parameters

    Will run the elbo plot to check quality across the parameters as well
    as save the lineage results.
    :param n: donor number
    :param indir:
    :param n_clone_list: clones parameter for vireo
    :param rerun_model: If file already there, whether to re-run or not
    :param prefix: Prefix for saving.
    :return:
    """
    ##TODO: change the import names in the old scripts
    cells_meta = join(indir, f"cells_meta.tsv")
    curr_ad_f = join(indir, f"donor{n}.AD.mtx")
    curr_dp_f = join(indir, f"donor{n}.DP.mtx")
    curr_labels_f = join(indir, f"donor{n}.labels.txt")
    curr_labels_rawID = join(indir, f"cell_labels.donor{n}.txt")
    print(curr_ad_f)
    print(curr_dp_f)
    curr_ad = mmread(curr_ad_f).tocsc()
    curr_dp = mmread(curr_dp_f).tocsc()
    curr_labels = pd.read_csv(curr_labels_f, index_col=0)
    curr_cells_meta = pd.read_csv(curr_labels_rawID, index_col=0)
    cells_meta = pd.read_csv(cells_meta, sep='\t')
    print(curr_cells_meta.head())
    run_elbo(curr_ad, curr_dp, join(outdir, f"donor{n}"), n_clone_list=n_clone_list,
             save_clusters=True, labels=curr_labels[["sample ID"]], rerun_model=rerun_model,
             cells_meta=curr_cells_meta)
    return


