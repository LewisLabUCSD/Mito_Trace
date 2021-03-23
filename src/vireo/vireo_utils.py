import pickle
from vireoSNP import Vireo
from vireoSNP.plot.base_plot import heat_matrix
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mplh.cluster_help as ch
import os
from os.path import join, exists

#prob_thresh = 0.9
#doublet_thresh = 0.9
#low_conf_cells = np.flatnonzero(doublet_prob > doublet_thresh)
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
        curr_sample_colors = sample_colors.iloc[cell_clusters[n]].copy()
        curr_sample_colors = curr_sample_colors.reset_index()

        cell_clusters[n] += 1
        # Get their IDs
        if outdir != "" and exists(outdir):
            curr_out = join(outdir, f"{out_f}_donor{n}_cells.txt")
            curr_str = "\n".join(cell_clusters[n].astype(str))
            with open(curr_out, "w") as f:
                f.write(curr_str)
            curr_sample_colors.to_csv(
                join(outdir, f"{out_f}_{n}.labels.txt"))
    return


### Run lineage tracing for each cluster individually
def run_elbo(ad, dp, out_f, sample_colors,
             n_clone_list=None, save_clusters=True,
             labels=None, rerun_model=False):
    if n_clone_list is None:
        n_clone_list = np.arange(2, 5)
    n_initials = 50
    _ELBO_mat = []

    for k in n_clone_list:
        print('lineages', k)
        curr_clone_model_f = out_f + f".clones{k}.modelCA.p"
        curr_lineage_f = out_f + f".clones{k}.lineages.png"
        if os.path.exists(curr_clone_model_f) and (not rerun_model):
            print('Already ran. Loading model')
            modelCA = pickle.load(open(curr_clone_model_f, "rb"))
        else:
            _models_all = []
            _elbo_temp = []
            for i in range(n_initials):
                _modelCA = Vireo(n_var=ad.todense().shape[0],
                                 n_cell=ad.todense().shape[1],
                                 n_donor=k, n_GT=2, fix_beta_sum=False,
                                 ASE_mode=True)
                _modelCA.set_prior(
                    beta_mu_prior=np.array([[0.01, 0.5]]))
                _modelCA.fit(ad, dp, min_iter=20, verbose=False)
                _elbo_temp.append(_modelCA.ELBO_[-1])
                _models_all.append(_modelCA)

            _ELBO_mat.append(_elbo_temp)
            _idx = np.argmax([x.ELBO_[-1] for x in _models_all])
            modelCA = _models_all[_idx]
            _losses = modelCA.ELBO_
            pickle.dump(modelCA,
                        open(out_f + f".clones{k}.modelCA.p", "wb"))


        # Run plot_vireo_out and extract the cell IDs as well
        ## Choose the model giving highest ELBO
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

    if len(_ELBO_mat) != 0:
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
        plt.show()
    return _ELBO_mat, n_clone_list


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

    fig = plt.figure(figsize=(7, 4), dpi=100)
    plt.subplot(1, 2, 1)
    im = heat_matrix(modelCA.ID_prob, cmap="Oranges", alpha=0.8,
                     display_value=False, row_sort=True)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.title("Assignment probability")
    plt.xlabel("Clone")
    plt.ylabel("%d cells" % (modelCA.n_cell))
    plt.xticks(range(modelCA.n_donor))
    plt.subplot(1, 2, 2)
    AF_SNPs = np.sum(
        modelCA.GT_prob * np.expand_dims(modelCA.beta_mu, 1), axis=2)
    if to_sqrt:
        AF_SNPs = np.sqrt(AF_SNPs)

    # f = plt.figure()
    im = heat_matrix(AF_SNPs, cmap="Blues", alpha=0.8,
                     display_value=False, row_sort=True)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.title("Mean allelic ratio")
    plt.xlabel("Clone")
    plt.ylabel("%d SNPs" % (modelCA.n_var))
    plt.xticks(range(modelCA.n_donor))
    plt.tight_layout()
    plt.savefig(out_f)
    plt.close()

    ch.plot_cluster(pd.DataFrame(AF_SNPs), cmap='Blues', alpha=0.8,
                    to_row_clust=True, to_col_clust=True,
                    to_legend=True, white_name=None)
    plt.savefig(out_f + ".variants.labels.png")
    # plt.show()
    return
