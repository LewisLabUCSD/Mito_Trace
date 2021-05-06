import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from os.path import join, basename, exists
import click
from glob import glob
from mplh.fig_utils import helper_save as hs
from mplh.cluster_help import plot_cluster as pc
import logging
import os
from matplotlib import rcParams
#import scanpy
#rcParams['lines.markersize'] = 10
import numpy as np


def load_lineage(in_files_dict, lin_col="Lineage"):
    """Load the cell indices for each cluster and each sample

    :param in_files_dict: Keys are names of cluster, values are the file paths
    :return dictionary where keys are names of cluster and values is list
    of cell names.
    """
    lineages = {}
    for i in in_files_dict:
        lineages[i] = pd.read_csv(in_files_dict[i])[["ID"]]
        lineages[i][lin_col] = i
    return pd.concat(lineages, ignore_index=True)


def load_clusters_and_tsne(indir):
    """Load Flt3l and control cluster assignments"""
    clusters = pd.read_csv(join(indir, "clustering", "graphclust", "clusters.csv"))
    tsne = pd.read_csv(join(indir, "tsne/2_components/projection.csv"))
    clusters["ID"] = clusters["Barcode"]
    tsne["ID"] = clusters["Barcode"]
    return clusters, tsne


def cell_labels_extract_clusters():
    # Combine the samples, create new cluster IDs (will merge another time)
    return


def cell_labels_extract_lineage():
    # Extract the cell names in each cluster, seeing which sample it comes
    # from
    return


def remap(clusters, samples):
    """
    :param clusters:
    :param lineages:
    :param samples: Ordered list of samples

    :return: clusters and lineages df but with ID updated to be same (lineage name is used).

    """

    clust_map = {i+1:samples[i] for i in range(len(samples))}
    clusters["Sample Num"] = clusters["ID"].apply(lambda x: x.split("-")[1])
    clusters["init ID"] = clusters["ID"].apply(
        lambda x: x.split("-")[0])
    clusters["merged ID"] = clusters["ID"].copy()

    # def change_cb(val):
    #     return val["init ID"] + "-1" + "_" + clust_map[val["Sample Num"]]
    clusters["ID"] = clusters.apply(lambda x: x["init ID"] + "-1_" + clust_map[int(x["Sample Num"])] , axis=1)
    # clusters.apply(lambda x: x["init ID"] + "-1_" + clust_map[
    #     int(x["Sample Num"])], axis=1)
    return clusters.set_index("ID")


def merge_all_dfs(clusters, lineages, tsne):
    """ Create df where each index is the cell, that merges all the results
    :return:
    """
    df = pd.merge(pd.merge(clusters, tsne, how='outer', left_index=True,
                      right_index=True), lineages, left_index=True,
             right_on="ID")
    df["Sample"] = df["ID"].apply(lambda x: x.split("_")[-1])
    return df


def count_overlap_lineage_cluster(df, out_f):
    """

    df: Columns: ["Cell", "Cluster", "Donor", "Lineage", "Sample]
    """
    df["ind"] = df["Lineage"].astype(str) + "_" + df["Donor"].astype(str)
    df["cluster"] = df["Cluster"].astype(str) + "_" + df["Sample"].astype(str)
    sns.countplot(y="ind", hue="cluster", data=df)
    hs(out_f)
    #plt.savefig(out_f, dpi=300)
    return


def plot_overlap_lineage_cluster(df, out_f, lin_col="Lineage"):
    print('overlap lineage', out_f)

    curr = df.groupby(["Cluster", lin_col]).size().reset_index().rename({0:"Count"}, axis=1)
    print(len(curr))
    if len(curr)==0:
        print('no matrix')
        return
    curr = curr.pivot(index="Cluster", columns=lin_col, values="Count").fillna(0).astype(int)
    sns.clustermap(curr, metric="cosine")
    hs(out_f)
    plt.close()
    return


def group_norm(df, anno_vals=("Cluster", "Sample")):
    anno_data = df[anno_vals].copy()
    df = df.drop(anno_vals,axis=1)/df.drop(anno_vals,axis=1).sum().sum()
    return pd.concat((anno_data, df),axis=1)


def group_transform(df, anno_vals=("Cluster", "Sample")):
    anno_data = df[anno_vals].copy()
    df = (df.drop(anno_vals,axis=1)+1)/(df.shape[0]+df.drop(anno_vals,axis=1).sum(axis=0)) # average over each lineage, not entire sample with a pseudocount
    return pd.concat((anno_data, df),axis=1)



def plot_percentage_cluster_enrich(cell_df, out_f=None, to_norm_allsample=True,
                                   all_donors=False, subplot=False):
    """ Similar to Fig. 2E in Dawn et al 2021 Nature Cell Bio

    Create the MTclone-by-FateCluster counts and run t-sne on this.
    Then Plots the t-sne and overlays the cluster counts on them.

    :param cell_df: pd.DataFrame where each row is cell and column is the cell attributes,
                    including cell type and lineage
    :return tsne_out
    :rtype scanpy object
    """

   # barcode_matrix = cell_df.groupby(["Cluster","Lineage"]).size().reset_index().rename({0:"Count"}, axis=1).pivot(index="Cluster", columns="Lineage", values="Count").fillna(0)
    if all_donors and to_norm_allsample: # if norm_allsample, will group by Donor as well to make sure normalized properly.
        anno_vals = ["Cluster", "Sample", "Donor"]
        norm_vals = ["Sample", "Donor"]
        #col_meta = cell_df.columns
    else:
        norm_vals = ["Sample"]
        anno_vals = ["Cluster", "Sample"]

    barcode_matrix = cell_df.groupby(anno_vals + ["Lineage"]).size().reset_index().rename(
        {0: "Count"}, axis=1).pivot_table(index=anno_vals,
                                          columns="Lineage",
                                          values="Count").fillna(0).reset_index()
    if all_donors:
        anno_coldata = pd.DataFrame(barcode_matrix.drop(anno_vals,axis=1).apply(lambda x: x.name.split('_')[1])).rename({0:"Donor"}, axis=1)
    else:
        anno_coldata = None

    if to_norm_allsample:
        barcode_matrix = barcode_matrix.groupby(
            norm_vals , group_keys=False).apply(group_norm, anno_vals)
    else:
        barcode_matrix = barcode_matrix.groupby(norm_vals,
            group_keys=False).apply(group_transform, anno_vals)
    # if 'Cluster' not in barcode_matrix.columns or "Sample" not in barcode_matrix.columns:
    #     print("Cluster and sample not in barcode_matrix")
    #     print('bc mat:')
    #     print(barcode_matrix.head())
    #     return
    barcode_matrix.to_csv(out_f+".csv")

    anno_data = barcode_matrix[anno_vals].copy()
    barcode_matrix = barcode_matrix.drop(anno_vals , axis=1)
    pc(barcode_matrix, row_meta=anno_data, to_col_clust=True, to_row_clust=True, metric="cosine", white_name=None, to_legend=True, col_meta=anno_coldata)
    hs(out_f+".png")
    plt.close()

    #if subplot:

    # Import into scanpy and run t-sne
    # import src.utils.scanpy_utils as scut
    # adata = scut.create_scanpy(df=barcode_matrix, sample_cols=barcode_matrix.columns)
    # scut.cluster_embedding(adata)
    # scut.umap_scanpy(adata, f_save=f_save)
    #create_scanpy(in_f=None,
    #              df: pd.DataFrame = None, sample_meta_f = None, sample_cols = None, is_header = True)
    return barcode_matrix


def plot_cluster_tsne(df, out_f, lin_col="Lineage", as_factor=False,
                      col_wrap=4):
    """
    Color by lineage
    :param df:
    :return:
    """
    print('tsne', out_f)

    if as_factor:
        print("here")
        g = sns.FacetGrid(df, col=lin_col,col_wrap=col_wrap)
        g.map_dataframe(sns.scatterplot, x="TSNE-1", y="TSNE-2", hue="Sample",
              s=20)
    else:
        sns.scatterplot(data=df.astype({lin_col:int}), x="TSNE-1", y="TSNE-2", hue =lin_col,
                        palette=sns.color_palette("Set2", len(set(df[lin_col]))), s=5)
                        #scatter_kws={"s": 10})
    hs(out_f, to_pdf=False, to_svg=False)
    plt.close()
    return


def plot_all_hex(df, col="Sample", plot_label=None,
                 out_f=None):
    #print('df',df.head())
    g = sns.FacetGrid(data=df, col=col)
    g.map(sns.kdeplot, "TSNE-1", "TSNE-2")
    for ax in g.axes.flat:
        # set aspect of all axis
        ax.set_aspect('equal', 'box')

    if plot_label is not None:
        plt.gcf().subplots_adjust(top=0.8)
        plt.suptitle(plot_label)

    if out_f != None:
        hs(out_f, to_svg=False, to_pdf=False)

    return df


def nonparametric_enrichment_overlap():
    return


def run_single_lineage(curr_lineage_files, tsne, clusters,  out_f, lin_col="Lineage",
                       to_plot=True, to_save=True):
    # 1. Get the lineage files and load them
    lineages = load_lineage(curr_lineage_files,
                            lin_col=lin_col)  # don_col= "Donor"

    # 2. Merge the dfs
    cell_df = merge_all_dfs(clusters, lineages, tsne)

    # 3. Save lineages in clusters
    # Save cell_df
    #out_f = join(outdir, f"cells_merged_donor{d}_nclones{n_clones}.csv")
    if to_save:
        cell_df.to_csv(out_f)
    # cell_df.to_csv(join(outdir, "cells_merged.csv"))

    if len(cell_df) == 0:
        logging.warning("No overlap here")
        return cell_df

    if to_plot:
        plot_percentage_cluster_enrich(cell_df, out_f=out_f+'.overlap'+'_percent_normAll',
                                       to_norm_allsample=True)
        plot_percentage_cluster_enrich(cell_df, out_f=out_f+'.overlap'+'_percent_normClone',
                                       to_norm_allsample=False)

        # 6. Plot the heatmap of overlap between lineage and cluster, across
        # the different donors
        plot_overlap_lineage_cluster(cell_df, out_f=out_f + ".overlap", lin_col=lin_col)

        # 7. Plot the tsne, but overlay the different lineges as different
        #    colors, across different donors.
        plot_cluster_tsne(cell_df, out_f + ".tsne.subplots",
                          lin_col=lin_col, as_factor=True, col_wrap=4)
        #plot_cluster_tsne(cell_df, out_f=out_f + ".tsne", lin_col=lin_col)

        plot_all_hex(cell_df, col="Sample", plot_label="T-SNE cell density",
                     out_f=out_f+".kde")

    return cell_df


@click.command()
@click.argument("aggregate_in", type=click.Path(exists=True))
@click.argument("lineage_in", type=click.Path(exists=True))
@click.argument("outdir", type=click.Path())
@click.argument("samples", type=click.STRING)
@click.argument("n_donors", type=click.INT)
@click.option("--only_donors", default=False)
@click.option("--nclones_values", default="")
def main(aggregate_in, lineage_in, outdir, samples, n_donors, only_donors=False,
         nclones_values=""):
    #aggregate_out
    logging.basicConfig()
    # A. Set up variables
    samples = samples.split(",")

    if nclones_values == "":
        nclones_values = []
    else:
        nclones_values = list(map(lambda x: int(x), nclones_values.split(",")))
    print('nclones_values', nclones_values)
    if not exists(outdir):
        import os
        print(f"making directory {outdir}")
        os.makedirs(outdir)

    # B. Load clusters and tsne, and map to samples
    clusters, tsne = load_clusters_and_tsne(aggregate_in)
    clusters = remap(clusters, samples=samples)
    tsne = remap(tsne, samples=samples)


    # If only_donors, just plot the multiplex results, not the subsequent
    # lineages.
    if only_donors:
        pair = dict()
        for d in range(n_donors):

            c = join(lineage_in,
                             f"cell_labels.donor{d}.txt")
            #curr_lineage_files.append(c)
            # c. add to the pair dictionary the {donor: file}
            pair[d] = c
        out_f = join(outdir,
                     f"cells_merged_lin_and_peak_donors.csv")
        run_single_lineage(pair, tsne, clusters, out_f, lin_col="Donor")
        return

    # C.
    # Construct the pairs dictionary which separates the files by
    # donors and n_clones param used in vireo.
    # The keys are (d, n_clones), and the value is a dictionary
    # where the keys are the lineage number and value is file name
    pair = dict()
    for d in range(n_donors):
        curr = glob(join(lineage_in, f"donor{d}_lineage_clones*.cell_labels.txt"))
        # break up by clones
        for c in curr:
            # a. get n_clones parameter used in vireo.
            curr_clones = int(basename(c).split("clones")[1].split("_")[0])

            if len(nclones_values) > 0 and curr_clones not in nclones_values:
                print(f"Not analyzing this nclone parameter: {curr_clones}")
                continue

            # b. get the cluster number
            curr_lin = basename(c).split("_")[-2].split(".")[0]

            # c. add to the pair dictionary
            if (d,curr_clones) not in pair:
                pair[(d, curr_clones)] = {curr_lin:c}
            else:
                pair[(d, curr_clones)][curr_lin] = c
            # if curr_clones not in clones:
            #     clones[curr_clones] = [c]
            # else:
            #     clones[curr_clones].append(c)

    # D.
    # For each donor-n_clones pair, merge the labels, plot results,
    # and save to csv.
    all_donors = {}
    for d, n_clones in pair:
        curr_lineage_files = pair[(d,n_clones)]
        if len(curr_lineage_files) < 3:
            logging.msg(f"Not enough files for donor {d} clones {n_clones}")
            continue
        out_f = join(outdir,
                     f"cells_merged_donor{d}_nclones{n_clones}.csv")
        curr_cell_df = run_single_lineage(curr_lineage_files, tsne, clusters, out_f,
                                          to_plot=False)
        # For each n_clones param, combine the donors
        curr_cell_df["Lineage"] = curr_cell_df["Lineage"].apply(lambda x: f"{x}_{d}")
        if n_clones in all_donors:
            all_donors[n_clones][d] = curr_cell_df
        else:
            all_donors[n_clones] = {d: curr_cell_df}

    for n_clone in all_donors:
        # subplots for each figure instead of being standalone, where
        # the dfs have an extra column Donor

        # Plot normalized lineage distribution with all donors
        all_donors_curr = pd.concat(all_donors[n_clone],sort=False).droplevel(1, axis=0).fillna(
            0).reset_index().rename({"index": "Donor"}, axis=1)

        out_f = join(outdir, f"cells_merged_lin_and_peak_nclones{n_clone}")
        plot_percentage_cluster_enrich(all_donors_curr, out_f=out_f+'.overlap'+'_percent_normClone',
                                    to_norm_allsample=False, all_donors=True)
        all_donors_curr.to_csv(out_f+"_allCells.csv")

    # pd.concat({"0": pd.DataFrame([[0, 1, 2], [4, 5, 6]],
    #                              columns=["A", "B", "C"]),
    #            "1": pd.DataFrame([[11, 0]], columns=["A", "C"])},
    #           sort=False).droplevel(1, axis=0).fillna(
    #     0).reset_index().rename({"index": "Donor"}, axis=1)

    return

# main(aggregate_in, lineage_in, outdir, samples, n_donors)


if __name__ == "__main__":
    lineage_in= "/data2/mito_lineage/Analysis/multiplex/data/jan21_2021/chrM/pseudo/minC200_minAF0.01/numC25000_ispropFalse/flt3/"
    aggregate_in  = "/data2/isshamie/mito_lineage/data/processed/mtscATAC/jan21_2021/MTblacklist/reanalysis/outs/analysis/"
    outdir = "/data2/mito_lineage/Analysis/lineage_and_peakclusters/results/jan21_2021/"
    samples = "J2,P2" #,J2"
    n_donors = '4'
    main([aggregate_in,
          lineage_in,
          outdir,
          samples,
          n_donors, "--only_donors", False,
          "--nclones_values", "20,100"])
    #    main(['aggregate_in', aggregate_in,
    #       'lineage_in',lineage_in,
    #       'outdir', outdir,
    #       'samples',samples,
    #       'n_donors',n_donors])
    #main()

