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


def filter_barcodes_from_CB(CB_read_number, cellr_bc_f):
    """
    Function that filters a dictionary where the keys are cells from a text file of cell barcodes, outputted from cellranger
    :param CB_read_number:
    :param cellr_bc_f:
    :return:
    """
    if cellr_bc_f is None:
        return CB_read_number
    filt_bc = pd.read_csv(cellr_bc_f, sep='\t', header=None)[0].values
    count = 0
    CB_read_number_filt = {}
    for x in filt_bc:
        if x not in set(CB_read_number.keys()):
            count += 1
            print(x)
        else:
            CB_read_number_filt[x] = CB_read_number[x]

    print(f"Number of missing barcodes in CB but present in filtered: {count}")
    print(f"Number of cells after using cellranger cell filter: {len(CB_read_number_filt)}")
    return CB_read_number_filt



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

