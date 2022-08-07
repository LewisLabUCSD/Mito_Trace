import scipy
from collections import defaultdict
import seaborn as sns
#mean_af_clust = mean_af.iloc[inds,cols]
import pandas as pd
import matplotlib.pyplot as plt
from dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage


def dendro_plot(df, row_meta):
    g = sns.clustermap(df, row_colors=row_meta,
                       row_cluster=True, col_cluster=True, vmax=0.2, vmin=0)
    inds = g.dendrogram_row.dendrogram["leaves"]
    #cols = g.dendrogram_col.dendrogram["leaves"]
    row_meta = row_meta.iloc[inds]
    return g, row_meta




def run_dynamic(df, metric='euclidean', method="average",
                minClusterSize=1, deepSplit=1,
                maxCoreScatter=None, minGap=None):
    if "multi" in df.columns.values:
        df = df.drop("multi", axis=1)
    distances = pdist(df.values, metric=metric)
    print('distances')
    print(distances)
    link = linkage(distances, method=method)
    print('link')
    print(link)
    print('deepsplit', deepSplit)
    clusters = cutreeHybrid(link, distances,
                            minClusterSize=minClusterSize,
                            deepSplit=deepSplit,
                            maxCoreScatter=maxCoreScatter,
                            minGap=minGap)["labels"]
    print('clusters', clusters)
    # clusters = pd.DataFrame({"ID":clusters["ID"], "labels": clusters["labels"]}, index=obj_norm_top10perc.index)#[["ID", "labels"]]
    clusters = pd.DataFrame({"labels": clusters},
                            index=df.index)  # [["ID", "labels"]]
    clusters["ID"] = df.index
    return clusters, link


def run_dynamic_hyper(df, metric='cosine', method="complete",
                      deepSplit=3):
    # clusters, link = run_dynamic(mean_af, metric='euclidean', method="average", minClusterSize=1)
    is_comp = False
    for i in ([3, 2, 1]):
        try:
            clusters, link = run_dynamic(df, metric='cosine',
                                         method="complete",
                                         minClusterSize=i, deepSplit=4)
            is_comp = True
            break
        except IndexError:
            print(f'didnt work for {i}')
    print("min cluster size:", i)

    if not is_comp:
        print('min size 1 and deepsplit 1')
        clusters, link = run_dynamic(df, metric=metric,
                                     method=method,
                                     minClusterSize=1, deepSplit=deepSplit)

    return clusters, link


def dendro_cluster(df, g=None, d_thresh=0.6, link=None, **kwargs):
    if g is not None:
        link = g.dendrogram_row.linkage
    if len(kwargs)>0:
        # Plotting function to truncate some if desired, but the output will be from the raw cell clusters

        f,ax = plt.subplots()
        scipy.cluster.hierarchy.dendrogram(link,
                                           labels = df.index, ax=ax,
                                           color_threshold=d_thresh, **kwargs)
    f2, ax2 = plt.subplots()
    den = scipy.cluster.hierarchy.dendrogram(link,
                                             labels = df.index, ax=ax2,
                                             color_threshold=d_thresh)
    plt.close(f2)
    return den



def old_v01_get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes


def old_v01_add_cluster_labels(df, den, row_meta):
    clusters = old_v01_add_cluster_labels(den)
    cluster = []
    for i in df.index:
        included=False
        for j in clusters.keys():
            if i in clusters[j]:
                cluster.append(j)
                included=True
        if not included:
            cluster.append("None")
    row_meta.loc[df.index,"den_clust"] = cluster
    return row_meta


def add_cluster_labels(den, row_meta):
    clusters = {}
    for ivl, lv_clust in zip(den["ivl"], den["leaves_color_list"]):
        #print(ivl, lv_clust)
        #print('ivl', ivl)
        #print('lv clr', lv_clust)
        clusters[ivl] = lv_clust
        #clusters[den['ivl'][lv]] = lv_clust

    out_meta = row_meta.copy()

    den_series = pd.Series(clusters)
    out_meta["den_clust"] = ""
    out_meta.loc[den_series.index, "den_clust"] = den_series

    return out_meta

def cluster_stats(clusters_df, out_f=None):
    from itertools import combinations
    size_pvals = []
    for pair in combinations(clusters_df["den_clust"].unique(), 2):
        print(pair)
        stat, p_val = scipy.stats.ranksums(clusters_df.loc[clusters_df["den_clust"] ==pair[0], "size"].values,
                                           y=clusters_df.loc[clusters_df["den_clust"] == pair[
                                                   1], "size"].values,
                                           alternative='two-sided')

        size_pvals.append(f"{pair[0]}, {pair[1]}, {str(p_val)}")

    if out_f is not None:
        with open(out_f + ".dendrogram_pvals.txt", "w") as f:
            # for l in size_pvals:
            f.write("\n".join(size_pvals))
    return size_pvals


def run_dendro_clust(df, dendroThresh, curr_clones, out_f=None,
                     to_clust_stats=False, use_seaborn=False, verbose=False, **kwargs ):
    #print('kwargs', kwargs)
    if use_seaborn:
        print('creating linkage with seaborn clustermap')
        g = sns.clustermap(data=df, row_cluster=True, col_cluster=True, **kwargs)
        inds = g.dendrogram_row.dendrogram["leaves"]
        cols = g.dendrogram_col.dendrogram["leaves"]
        curr_clones = curr_clones.iloc[inds]
        plt.close(g.fig)
        den = dendro_cluster(df, g, d_thresh=dendroThresh)
    else:
        print('creating linkage with scipy.cluster.hierarchy.linkage')
        from scipy.cluster.hierarchy import dendrogram, linkage
        from matplotlib import pyplot as plt
        linked = linkage(df, 'average')
        plt.figure(figsize=(10, 7))
        den = dendro_cluster(df, g=None, link=linked, d_thresh=dendroThresh)
        # den = dendrogram(linked, orientation='top', labels=df.index,
        #            distance_sort='descending', show_leaf_counts=True)

    if verbose:
        print('icoord and color list')
        print(len(den["icoord"]))
        print(len(den["color_list"]))
        print('leaves leaves color and df shape')
        print(len(den["leaves"]))
        print(len(den["leaves_color_list"]))
        print(df.shape)


    # cluster_classes = dc.get_cluster_classes(den)
    curr_clones = add_cluster_labels(den, curr_clones)

    if to_clust_stats:
        size_pvals = cluster_stats(curr_clones, out_f= out_f )

    #curr_clones["donor"]= [x.split("_")[0] for x in curr_clones.index]
    curr_clones = curr_clones.loc[df.index]

    # with open(out_f + ".dendrogram_pvals.txt", "w") as f:
    #     #for l in size_pvals:
    #     f.write("\n".join(size_pvals))
    # import mplh.cluster_help as ch
    # g = ch.plot_cluster(df, row_meta=curr_clones,
    #                     to_row_clust=True, to_col_clust=True,
    #                     row_clr_schemes={"size": "sequential",
    #                                      "donor": "categorical",
    #                                      "den_clust": "categorical"})
    if to_clust_stats:
        return den, curr_clones, size_pvals
    return den, curr_clones


