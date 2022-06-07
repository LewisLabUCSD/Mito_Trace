import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.style.use('fivethirtyeight')
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import scipy.stats as stats
import numpy as np
from scipy.stats import zscore
mpl.use('Agg')
from mplh.fig_utils import helper_save
import tqdm
from os.path import join
from mplh.fig_utils import helper_save as hs
from pandarallel import pandarallel




def get_continuous_colors(df, col, clr_key=1, clr_type='sequential'):
    anno_labels = np.sort(df[col].unique())


    clr_keys = {1:sns.cubehelix_palette(len(anno_labels),
                                        light=.9, dark=.2, reverse=True,
                                        rot=.1, start=2.8),
                2:sns.cubehelix_palette(len(anno_labels),
                                        light=.9, dark=.2, reverse=True,
                                        rot=.1, start=4.2),
                3:sns.cubehelix_palette(len(anno_labels),
                        light=.9, dark=.2, reverse=True,
                        rot=.3, start=0)}
    divergent_clr_keys = {1:sns.diverging_palette(240, 10, n=len(anno_labels)),
                          2:sns.diverging_palette(150, 275, s=80, l=55, n=len(anno_labels))
                          }
    cat_clr_keys = {1:sns.color_palette("Set2"), 
                    2:sns.color_palette("Paired")
                    }
    if clr_type=="sequential":
        anno_pal = clr_keys[clr_key]
    elif clr_type=="divergent":
        anno_pal = divergent_clr_keys[clr_key]
    elif clr_type=="categorical":
        anno_pal = cat_clr_keys[clr_key]
    else:
        raise ValueError

    anno_lut = dict(zip(map(str, anno_labels), anno_pal))

    anno_colors = pd.Series(anno_lut)
    anno_colors

    df[f"{col}_map"] = df[col].apply(lambda x: anno_colors.loc[str(x)])

    return df, anno_labels, anno_lut

def plot_continuous_legend(g, anno_labels, anno_lut, n_labs=-1, title=None, loc='right', decimal_places=3):
    if n_labs == -1 or n_labs>len(anno_labels):
        n_labs = len(anno_labels)
    
    step = int(np.round(len(anno_labels)/n_labs))

    # if title is not None:
    #     tit = g.ax_heatmap.bar(0, 0,label=title, linewidth=0)
    for label in anno_labels[::step]:   
        if type(label) == str:
            g.ax_heatmap.bar(0, 0, color=anno_lut[str(label)],
                                    label=f'{label}', linewidth=0)  
        else:            
            g.ax_heatmap.bar(0, 0, color=anno_lut[str(label)],
                                    label=f'{label:.3g}', linewidth=0)
        # g.ax_col_dendrogram.bar(0, 0, color=anno_lut[str(label)],
        #                         label=label, linewidth=0)
        # plt.bar(0, 0, color=anno_lut[str(label)],
        #                 label=label, linewidth=0)

    g.ax_heatmap.legend(bbox_to_anchor=(1.4, 1.2), ncol=4, loc=loc, borderaxespad=1)
    #g.ax_col_dendrogram.legend(ncol=4 )#loc="best", ncol=6)
    return g.ax_heatmap.legend()



def biplot(score, coeff, labels=None, PCs=(1, 2)):
    xs = score[:, 0]
    ys = score[:, 1]
    n = coeff.shape[0]
    scalex = 1.0 / (xs.max() - xs.min())
    scaley = 1.0 / (ys.max() - ys.min())
    plt.scatter(xs * scalex, ys * scaley)
    for i in range(n):
        plt.arrow(0, 0, coeff[i, 0], coeff[i, 1], color='r',
                  alpha=0.5)
        if labels is None:
            plt.text(coeff[i, 0] * 1.15, coeff[i, 1] * 1.15,
                     "Var" + str(i + 1), color='g', ha='center',
                     va='center')
        else:
            plt.text(coeff[i, 0] * 1.15, coeff[i, 1] * 1.15,
                     labels[i], color='g', ha='center', va='center')
            plt.xlim(-1, 1)
            plt.ylim(-1, 1)
            plt.xlabel("PC{}".format(PCs[0]))
            plt.ylabel("PC{}".format(PCs[1]))
            plt.grid()


def construct_pca(data, top_comp=5, num_feat_per=3, num_comp=4, save_f=None):
    p = PCA()
    p.fit(data)
    embedding = p.transform(data)
    variance = p.explained_variance_ratio_  # calculate variance ratios
    plt.figure()
    var = np.cumsum(
        np.round(p.explained_variance_ratio_, decimals=3) * 100)
    plt.ylabel('% Variance Explained')
    plt.xlabel('# of Features')
    plt.title('PCA Analysis')
    plt.ylim(30, 100.5)
    plt.scatter(x=np.arange(1, len(var) + 1), y=var)
    if save_f is not None:
        helper_save(save_f + '_pcVariance.png')

    data["embedding_1"] = embedding[:, 0]
    data["embedding_2"] = embedding[:, 1]
    # Call the function. Use only the 2 PCs.

    plt.figure()
    biplot(data[["embedding_1", "embedding_2"]].values,
           np.transpose(p.components_[0:2, :]),
           labels=data.columns.values)
    
    if save_f is not None:
        helper_save(save_f + '_biplot.png')

    data = data.drop(["embedding_1", "embedding_2"], axis=1)
    inds_to_keep = set()
    for i in range(top_comp):
        # Get the top num_feat_per variants of the i'th component
        inds_to_keep = inds_to_keep.union(
            set((p.components_[i, :]).argsort()[::-1][:num_feat_per]))

    inds_to_keep = list(inds_to_keep)

    data = data.iloc[:, inds_to_keep]

    p_comp = p.components_[:num_comp, inds_to_keep]

    # Sort by top component
    top_c = p_comp[0, :].argsort()[::-1]
    p_comp = p_comp[:, top_c]
    data = data.iloc[:, top_c]


    #Transpose for plotting
    p_comp = p_comp.transpose()

    f, ax = plt.subplots(nrows=int(np.ceil(p_comp.shape[1] / 2)),
                         ncols=2, squeeze=True, figsize=(10, 10),
                         sharex='all') #, sharey='all'

    for feat in range(p_comp.shape[1]):
        x = np.arange(p_comp.shape[0])
        xl = data.columns.values
        if feat == 0:

            ax[0, 0].bar(x=x,
                         height=p_comp[:, feat], )
            #ax[0, 0].set_title(data.columns.values[feat])
            ax[0, 0].set_title(f'Component #{feat}')

            ax[0, 0].set_xticks(x)
            ax[0, 0].set_xticklabels(xl, rotation=90)

            ax[0, 0].set_xlabel("Variant")
            ax[0, 0].set_ylabel("Feature Load")
        elif int(feat / 2) == ax.shape[0]:
            continue
        else:

            ax[int(feat / 2), feat % 2].bar(
                x=x, height=p_comp[:, feat])
            # ax[int(feat / 2), feat % 2].set_title(
            #     data.columns.values[feat])
            ax[int(feat / 2), feat % 2].set_title(f'Component #{feat}')

            ax[int(feat / 2), feat % 2].set_xlabel("Variant")
            ax[int(feat / 2), feat % 2].set_ylabel(
                "Feature Load")  # ax[int(feat / 2), feat % 2].set_xticklabels()

            ax[int(feat / 2), feat % 2].set_xticks(x)
            ax[int(feat / 2), feat % 2].set_xticklabels(xl, rotation=90)


    plt.subplots_adjust(hspace=0.4)        # plt.legend()
    #helper_save("PCloadings")
    if save_f is not None:
        helper_save(save_f + '_pcLoadings.png')
    return data, p_comp


def expand_df(df, results_expand, results_key, expand_col):
    """ Expanding a df by duplicating each row and adding values from results_expand

    :param df: DataFrame to expand. Each index needs to be a key in results_expans
    :param results_expand: The dictionary where the results are going to be used for each row
    :param results_key: The key in results_expand that have the values to add
    :param expand_col: The column to add that will contain the list values
    :return:
    """
    df_full = []
    for ind, val in df.iterrows():
        curr_results = results_expand[ind]
        b_a_df = curr_results[results_key]
        # Each element will be a list of the values.
        s = df.loc[ind].copy()
        curr = pd.DataFrame([s.copy()] * len(b_a_df))
        curr[expand_col] = b_a_df[expand_col].values
        df_full.append(curr)

    df_full = pd.concat(df_full)
    return df_full





def compare_arbitrary_labels(l1, l2):
    """ Bring cluster labels into true labels by looking into overlap between them

    Two lists of categorical labels. It will try to put together the same
    labels for each using brute force, by first calling the cluster with
    the most matches for a label, and then go one by one till all labels
    are used. Any additional ones will be added ints
    :param l1: Each element is a label for that index.
    :type l1: iterable
    :param l2: Same length of l1, and each element is a label for that index.
    :return:
    :rtype list
    """
    l1 = np.array(l1)
    l2 = np.array(l2)
    # Get labels from l1 and l2
    l1_labels = list(set(l1))
    l2_labels = list(set(l2))

    # For each l1 label, create a map of overlap to the l2 labels
    df = pd.DataFrame(index=l1_labels, columns=l2_labels, dtype=int)
    df.loc[:,:] = 0
    for ind, curr_l1 in enumerate(l1):
        curr_l2 = l2[ind]
        df.at[curr_l1, curr_l2] += 1

    # Create label map
    l2_to_l1 = dict()
    while len(df.columns) > 0 and len(df.index)>0:
        rows, cols = np.where(df.to_numpy() == np.max(df.to_numpy()))
        # If many rows and cols, take the first
        l2_to_l1[df.columns[cols[0]]] = df.index[rows[0]]
        df = df.drop(df.index[rows[0]], axis=0).drop(df.columns[cols[0]], axis=1)
    # If there are still labels in l2, make them labelled as ints based on the number of labels
    if len(df.columns) > 0 and len(df.index) == 0:
        for c in df.columns:
            curr = len(l2_to_l1) + 1
            while curr in l2_to_l1.values(): # in case it was already labelled
                curr += 1
            l2_to_l1[c] = curr
    print(l2_to_l1)
    l2_new = list(map(lambda x: l2_to_l1[x], l2))
    return l2_new


def get_frag(indir):
    df = pd.read_csv(join(indir, "fragments.tsv.gz"), sep="\t")
    return


def check_frag_in_peak(frag, peaks):
    # print(((frag["Start"] <= peaks["End"]) & (frag["End"] >= peaks["Start"])).any())
    return ((frag["Start"] <= peaks["End"]) & (
                frag["End"] >= peaks["Start"])).any()


def calc_frip_per_cell(frags, peaks=None, cells=None, to_annotate=False,
                       to_parallel=False, nb_workers=12):
    frags = (frags.copy())[frags["Cell"].isin(cells)]
    if to_annotate:
        if to_parallel:
            pandarallel.initialize(nb_workers=nb_workers)
            frags["in_peak"] = frags.parallel_apply(check_frag_in_peak,
                                                    args=(peaks,),
                                                    axis=1)
        else:
            frags["in_peak"] = frags.apply(check_frag_in_peak, args=(peaks,), axis=1)
    cb_frags = frags.groupby("Cell")
    frags_sum = cb_frags["Count"].sum()
    frags_in_read_sum = \
    frags.groupby("in_peak").get_group(True).groupby("Cell")[
        "Count"].sum()
    frip = frags_in_read_sum / frags_sum
    frip_df = pd.concat((frags_sum.rename("Total frags"),
                         frags_in_read_sum.rename("Frags in peak"),
                         frip.rename("FRIP")), axis=1)
    return frip_df


def plot_frip_per_cell(frip_df, out_f=None):
    f, ax = plt.subplots(nrows=2, ncols=1)
    sns.violinplot(frip_df, y="frip", ax=ax[0])
    sns.scatterplot(frip_df, x="frags_sum", y="frip", ax=ax[1])
    hs(out_f)
    return



def test_compare_arbitrary_labels():
    l1 = [0, 1, 2, 2, 2, 2, 1, 0, 0, 3, 0]
    l2 = [2, 0, 1, 2, 8, 8, 0, 7, 2, 7, 2]
    l2_new = compare_arbitrary_labels(l1,l2)
    print(l2_new)

    l1 = [2, 0, 1, 2, 8, 8, 0, 7, 2, 7, 2]
    l2 = [0, 1, 2, 2, 2, 2, 1, 0, 0, 3, 0]
    l2_new = compare_arbitrary_labels(l1,l2)
    print(l2_new)



if __name__ == '__main__':
    test_compare_arbitrary_labels()
