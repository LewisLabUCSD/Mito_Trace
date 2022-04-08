import scipy
from collections import defaultdict
import seaborn as sns
#mean_af_clust = mean_af.iloc[inds,cols]


def dendro_plot(df, row_meta):
    g = sns.clustermap(df, row_colors=row_meta,
                       row_cluster=True, col_cluster=True, vmax=0.2, vmin=0)
    inds = g.dendrogram_row.dendrogram["leaves"]
    #cols = g.dendrogram_col.dendrogram["leaves"]
    row_meta = row_meta.iloc[inds]
    return g, row_meta


def dendro_cluster(df, g, d_thresh=0.6):
    den = scipy.cluster.hierarchy.dendrogram(g.dendrogram_row.linkage,
                                             labels = df.index,
                                             color_threshold=d_thresh)
    return den


def get_cluster_classes(den, label='ivl'):
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


def add_cluster_labels(df, den, row_meta):
    clusters = get_cluster_classes(den)
    cluster = []
    for i in df.index:
        included=False
        for j in clusters.keys():
            if i in clusters[j]:
                cluster.append(j)
                included=True
        if not included:
            cluster.append(None)
    row_meta.loc[df.index,"den_clust"] = cluster
    return row_meta


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
        with open(out_f + "dendrogram_pvals.txt", "w") as f:
            # for l in size_pvals:
            f.write("\n".join(size_pvals))
    return size_pvals
