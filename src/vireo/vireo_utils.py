import pickle
from vireoSNP import Vireo
from vireoSNP.plot.base_plot import heat_matrix
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mplh.cluster_help as ch
import os
from os.path import join, exists
from joblib import Parallel, delayed
from scipy.io import mmread
from src.pseudo_batch import wrap_write_mtx_df


def run_single(ad_dp, n):
    _modelCA = Vireo(n_var=ad_dp[0].todense().shape[0],
                     n_cell=ad_dp[0].todense().shape[1], n_donor=n, n_GT=2,
                     fix_beta_sum=False, ASE_mode=True)
    _modelCA.set_prior(beta_mu_prior=np.array([[0.01, 0.5]]))
    _modelCA.fit(ad_dp[0], ad_dp[1], min_iter=20, verbose=False)
    return _modelCA


def run_vireo(ad, dp, k, out_f=None, plot_qc=False, n_cores=24,
              n_initials=50):
    _ELBO_mat = []
    # _models_all = []
    # _elbo_temp = []
    # for i in range(n_initials):
    #     print(i)
    #     _modelCA = Vireo(n_var=ad.todense().shape[0],
    #                      n_cell=ad.todense().shape[1],
    #                      n_donor=k, n_GT=2, fix_beta_sum=False,
    #                      ASE_mode=True)
    #     _modelCA.set_prior(
    #         beta_mu_prior=np.array([[0.01, 0.5]]))
    #     _modelCA.fit(ad, dp, min_iter=20, verbose=False)
    #     _elbo_temp.append(_modelCA.ELBO_[-1])
    #     _models_all.append(_modelCA)
    _models_all = Parallel(n_jobs=n_cores)(delayed(run_single)((ad, dp), k) for i in range(n_initials))

    _elbo_temp = [x.ELBO_[-1] for x in _models_all]
    _idx = np.argmax([x.ELBO_[-1] for x in _models_all])
    modelCA = _models_all[_idx]
    _losses = modelCA.ELBO_

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
                     sample_colors,
                     outdir=None, out_f="lineage"):
    low_conf_cells = np.flatnonzero(doublet_prob > doublet_thresh)
    cell_clusters = dict()
    # cell_clusters_names = dict()
    for n in range(modelCA.ID_prob.shape[1]):
        # Drop low probability and/or high doublet probability
        cell_clusters[n] = np.flatnonzero(
            (modelCA.ID_prob[:, n] > prob_thresh))
        cell_clusters[n] = cell_clusters[n][
            ~(np.isin(cell_clusters[n], low_conf_cells))]

        if sample_colors is not None:
            curr_sample_colors = sample_colors.iloc[cell_clusters[n]].copy()
            curr_sample_colors = curr_sample_colors.reset_index()

        cell_clusters[n] += 1
        # Get their IDs
        if outdir != "" and exists(outdir):
            curr_out = join(outdir, f"{out_f}_donor{n}_cells.txt")
            curr_str = "\n".join(cell_clusters[n].astype(str))
            with open(curr_out, "w") as f:
                f.write(curr_str)
            if sample_colors is not None:
                curr_sample_colors.to_csv(
                    join(outdir, f"{out_f}_{n}.labels.txt"))
    return cell_clusters


def separate_donors(AD, DP, modelCA, sample_labels, OUTDIR, N_DONORS,
                    doublet_prob, sample_colors,
                    prob_thresh=0.9, doublet_thresh=0.9):
    """Separates the matrices and labels by donor using the multiplex output.

    :param AD:
    :param DP:
    :param modelCA:
    :param sample_labels:
    :param OUTDIR:
    :param N_DONORS:
    :param doublet_prob:
    :param sample_colors:
    :param prob_thresh:
    :param doublet_thresh:
    :return:
    """
    low_conf_cells = np.flatnonzero(doublet_prob > doublet_thresh)
    cell_clusters = dict()

    # For each donor, extract their cells for their matrices and labels
    for n in range(N_DONORS):
        # Drop low probability and/or high doublet probability
        cell_clusters[n] = np.flatnonzero(
            (modelCA.ID_prob[:, n] > prob_thresh))
        cell_clusters[n] = cell_clusters[n][
            ~(np.isin(cell_clusters[n], low_conf_cells))]

        # Get cell labels
        curr_sample_colors = sample_colors.iloc[cell_clusters[n]].copy()
        curr_sample_colors = curr_sample_colors.reset_index()
        curr_sample_labels = sample_labels.iloc[cell_clusters[n]]

        # Change the index to the sparse matrix index
        curr_ad = pd.DataFrame(
            AD.todense()[:, cell_clusters[n]]).reset_index().melt(
            id_vars='index', var_name="Cell",
            value_name="Count").rename({"index": "Position"}, axis=1)
        curr_dp = pd.DataFrame(
            DP.todense()[:, cell_clusters[n]]).reset_index().melt(
            id_vars='index', var_name="Cell",
            value_name="Count").rename({"index": "Position"}, axis=1)
        # Drop 0s to make sparse
        curr_ad = curr_ad.loc[~(curr_ad["Count"] == 0)]
        curr_dp = curr_dp.loc[~(curr_dp["Count"] == 0)]

        # Update the cell and position maps, making it 1-based
        curr_cell_map = {val: ind + 1 for ind, val in
                         enumerate(np.sort(curr_dp["Cell"].unique()))}
        curr_pos_map = {val: ind + 1 for ind, val in enumerate(
            np.sort(curr_dp["Position"].unique()))}
        # curr_ad["Cell"].map(curr_cell_map)
        curr_ad["Cell"] = curr_ad["Cell"].map(curr_cell_map)
        curr_ad["Position"] = curr_ad["Position"].map(curr_pos_map)
        curr_dp["Cell"] = curr_dp["Cell"].map(curr_cell_map)
        curr_dp["Position"] = curr_dp["Position"].map(curr_pos_map)

        print(f"Donor {n}: {len(cell_clusters[n])} cells ")
        print(curr_dp.shape)
        print(curr_ad.shape)

        cell_clusters[n] += 1
        # Get their IDs
        if OUTDIR != "" and exists(OUTDIR):
            curr_out = join(OUTDIR, f"donor{n}_cells.txt")
            curr_str = "\n".join(cell_clusters[n].astype(str))
            with open(curr_out, "w") as f:
                f.write(curr_str)
            curr_sample_colors.to_csv(
                join(OUTDIR, f"donor{n}.labels.txt"))
            curr_sample_labels.to_csv(
                join(OUTDIR, f"cell_labels.donor{n}.txt"))
            print('curr_ad', curr_ad.head())
            wrap_write_mtx_df(OUTDIR, curr_ad, curr_dp,
                              oth=None, to_rm=True, prefix=f"donor{n}",
                              columns=('Position', 'Cell', 'Count'))
    if OUTDIR != "":
        AF_SNPs = np.sum(
            modelCA.GT_prob * np.expand_dims(modelCA.beta_mu, 1),
            axis=2)
        pd.DataFrame(AF_SNPs, columns=[f"Cluster {x}" for x in
                                       range(AF_SNPs.shape[1])]).to_csv(
            join(OUTDIR, "AF_SNPs.csv"), index=False)
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
        ch.plot_cluster(clust_df, cmap='Oranges', alpha=0.8,
                        to_row_clust=True, to_col_clust=False,
                        row_meta=labels, to_legend=True,
                        white_name=None)
        plt.suptitle("Cell-cluster probability")
        plt.savefig(out_f + ".labels.png")  # plt.close()


    #f, ax = plt.subplots(figsize=(7, 4), dpi=300, nrows=2, ncols=2)
    # plt.figure()
    # plt.subplot(1, 2, 1)
    # im = heat_matrix(modelCA.ID_prob, cmap="Oranges", alpha=0.8,
    #                  display_value=False, row_sort=True)
    # plt.colorbar(im, fraction=0.046, pad=0.04)
    # plt.title("Assignment probability")
    # plt.xlabel("Clone")
    # plt.ylabel("%d cells" % (modelCA.n_cell))
    # plt.xticks(range(modelCA.n_donor))
    #plt.subplot(1, 2, 2)
    AF_SNPs = np.sum(
        modelCA.GT_prob * np.expand_dims(modelCA.beta_mu, 1), axis=2)

    ch.plot_cluster(pd.DataFrame(AF_SNPs), cmap='Blues', alpha=0.8,
                    to_row_clust=True, to_col_clust=True,
                    to_legend=True, white_name=None)
    plt.savefig(out_f + ".variants.labels.png")
    # plt.show()
    if doublet_prob is not None and labels is not None:
        return clust_df, AF_SNPs
    return None


### Run lineage tracing for each cluster individually
def run_elbo(ad, dp, out_f, sample_colors,
             n_clone_list=None, save_clusters=True,
             labels=None, rerun_model=True, n_cores=16, n_initials=50):
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
            doublet_prob = \
            modelCA.predict_doublet(ad, dp, update_GT=False,
                                    update_ID=False)[0].sum(axis=1)
            plot_vireo_out(modelCA, out_f=curr_lineage_f, to_sqrt=True,
                           labels=labels, doublet_prob=doublet_prob)

            extract_clusters(modelCA, prob_thresh=0.9,
                             doublet_thresh=0.9, doublet_prob=doublet_prob,
                             sample_colors=sample_colors,
                             outdir=out_f, out_f="lineage")

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
    curr_ad_f = join(indir, f"donor{n}.AD.mtx")
    curr_dp_f = join(indir, f"donor{n}.DP.mtx")
    curr_labels_f = join(indir, f"donor{n}.labels.txt")
    curr_labels_rawID = join(indir, f"cell_labels.donor{n}.txt")
    print(curr_ad_f)
    print(curr_dp_f)
    curr_ad = mmread(curr_ad_f).tocsc()
    curr_dp = mmread(curr_dp_f).tocsc()
    curr_labels = pd.read_csv(curr_labels_f, index_col=0)
    curr_sample_labels = pd.read_csv(curr_labels_rawID, index_col=0)
    print(curr_sample_labels.head())
    run_elbo(curr_ad, curr_dp, join(outdir, f"donor{n}"), n_clone_list=n_clone_list,
             save_clusters=True, labels=curr_labels[["sample ID"]], rerun_model=rerun_model,
             sample_colors=curr_sample_labels)
    return

