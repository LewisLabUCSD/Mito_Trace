import scanpy as sc
import pandas as pd
import anndata as ad


def read_meta(adata, in_f=None, df=None, cols=None, is_header=True,
              samp_or_feat='sample'):
    if not (in_f is None):
        if is_header:
            anno = pd.read_csv(in_f, sep='\t')
        else:
            anno = pd.read_csv(in_f, sep='\t', header=None)
    elif not (df is None):
            anno = df
    else:
        print("No meta")
        return adata
    if cols is not None:
        for c in cols:
            if samp_or_feat == "sample":
                adata.obs[c] = anno[c]
            if samp_or_feat == "feature":
                adata.var[c] = anno[c]
    else:
        for c in anno.columns:
            if samp_or_feat == "sample":
                adata.obs[c] = anno[c]
            elif samp_or_feat == "feature":
                adata.var[c] = anno[c]
    return adata


def create_scanpy(in_f=None, df:pd.DataFrame=None,
                  sample_meta_f=None, sample_df=None,
                  sample_cols=None, sample_is_header=True,
                  feature_meta_f=None, feature_df=None,
                  feature_cols=None, feature_is_header=True):
    """ Creates a scanpy object, either by passing in a filename or a
    pd.DataFrame. Additionally takes in sample meta-information based
    on either sample_cols which are columns from the dataframe,
    or from a separate csv file. Can take in feature information as well.

    :param in_f:
    :param df:
    :param sample_meta_f:
    :param sample_cols:
    :param is_header:
    :return:
    """
    if not (in_f is None):
        adata = sc.read(in_f)
    elif df is not None:
        adata = ad.AnnData(df) #.to_numpy(), obs=df.index, var=df.columns, dtype='int32')
    else:
        print('no file or df found')
        return

    ## Sample information
    adata = read_meta(adata, in_f=sample_meta_f, df=sample_df,
                      cols=sample_cols, is_header=sample_is_header,
                      samp_or_feat='sample')

    adata = read_meta(adata, in_f=feature_meta_f, df=feature_df,
                      cols=feature_cols, is_header=feature_is_header,
                      samp_or_feat='feature')
    return adata


def cluster_embedding(adata):
    sc.tl.leiden(adata)


def umap_scanpy(adata, color=None, f_save=None):
    sc.tl.umap(adata)
    if color is None:
        sc.pl.umap(adata, color=adata.var_names)
    else:
        sc.pl.umap(adata, color=color)
    if f_save is not None:
        adata.write(f_save)
        #import matplotlib.pyplot as plt
        #plt.savefig(f_save)
    return adata


def test_create_scanpy():
    print('here')
    af_f = "/data2/mito_lineage/data/processed/mttrace/2020_11_18/PBMC_J/mapq_0/cellr_True/PBMC_J_200/PBMC_J.af.tsv"
    sample_meta_f = "/data2/mito_lineage/data/processed/mttrace/2020_11_18/PBMC_J/mapq_0/cellr_True/PBMC_J_200/PBMC_J.depthTable.txt"
    adata = create_scanpy(af_f, sample_meta_f=sample_meta_f,
                          sample_cols=[0], is_header=False)
    return adata


print('here')
if __name__ == "__main__":
    test_create_scanpy()

