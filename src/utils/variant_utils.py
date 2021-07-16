import pandas as pd
import numpy as np
import os


def preprocess_variants(variants, style=">"):
    if style == ">":
        def split(x):
            s = x.split(">")
            return [s[0][:-1], s[0][-1], x[-1]]

        curr = pd.DataFrame(list(map(split, variants)),
                            columns=["position", "ref", "alt"],
                            index=variants)
    else:
        print(f"style {style} not implemented")
        return
    return curr


def annotate_vcf(vcf_in, gff_in, out_f):
    cmd = f"bedtools annotate -i {vcf_in} -files {gff_in} > {out_f}"
    print(cmd)
    os.system(cmd)
    return


def variants_dense(AF_df, vars_to_plot, samples_d, donors_d,
                   variant_d=None, lineage_d=None):
    #var_sort = AF_df.mean(axis=1).argsort()[::-1]
    AF_df = AF_df.rename_axis("Cell", axis=1).rename_axis("Variant", axis=0)
    if vars_to_plot is None:
        vars_to_plot=len(AF_df)
    var_sort = AF_df.mean(axis=1).sort_values()[::-1].index[:vars_to_plot]
    variants_box=AF_df.loc[var_sort].reset_index().melt(id_vars='Variant',
                                                                              value_name='AF')
    variants_box['condition'] = variants_box['Cell'].map(samples_d)
    variants_box['donor'] = variants_box['Cell'].map(donors_d)
    #variants_box=variants_box.dropna().copy()
    variants_box['sqrtAF'] = np.sqrt(variants_box['AF'])
    if variant_d is not None:
        variants_box['variant type'] = variants_box['Variant'].map(variant_d)
    if lineage_d is not None:
        variants_box['lineage'] = variants_box['Cell'].map(lineage_d)
    return variants_box


def type_of_variants(variants, style=">", to_preproc=True):
    if to_preproc:
        variants = preprocess_variants(variants, style=style)
    # Get types of mutations
    def var_type(x):
        nts = set(x[["ref", "alt"]])
        if "N" in nts:
            return "Undefined"
        if nts == {"A","G"} or nts == {"T", "C"}:
            return "Transition"
        return "Transversion"
    variants["variant type"] = variants.apply(var_type, axis=1)

    variants["variant change"] = variants["ref"]+">"+variants["alt"]
    return variants


def filt_high(df, thresh):
    """ Remove rows that have an average higher than threshold

    :param df:
    :param thresh:
    :return:
    """
    return df.loc[(df.mean(axis=1)<thresh)].index.values


def filt_low(df, thresh):
    """ Remove rows that have an average lower than threshold

    :param df:
    :param thresh:
    :return:
    """
    return df.loc[(df.mean(axis=1)>thresh)].index.values


def extract_variants_per_donor(af, cell_meta, af_thresh=0, cells_thresh=1):
    af = af>af_thresh
    samps = set(cell_meta.values)
    variants_dict = {}
    for s in samps:
        print(af.loc[(af.loc[:,cell_meta[cell_meta==s].index.astype(object)]>af_thresh).any(axis=1)].index)
        variants_dict[s] =  af.loc[(af.loc[:,cell_meta[cell_meta==s].index.astype(object)]>af_thresh).any(axis=1)].index
    return variants_dict
