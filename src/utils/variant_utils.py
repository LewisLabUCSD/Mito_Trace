import pandas as pd
import numpy as np
import os
import seaborn as sns


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


def calc_variants_qc(other_af, other_dp, prefix=""):
    thresh_vals = [0.01, 0.05, 0.1, 0.9]
    thresh_cols = ["perc", "number"]
    qc_d = pd.DataFrame(index=other_af.index,
                        columns=["group", "lineage", "mean", "std", "depth mean",
                                 "depth number_1", "depth number_5",
                                 "depth number_10"] + [f"{pref}_{t}" for
                                                       t in thresh_vals
                                                       for pref in
                                                       thresh_cols])
    for t in thresh_vals:
        qc_d[f"perc_{t}"] = (other_af > t).sum(axis=1) / other_af.shape[1]
        qc_d[f"number_{t}"] = (other_af > t).sum(axis=1)
    qc_d["mean"] = other_af.mean(axis=1)
    qc_d["std"] = other_af.std(axis=1)
    qc_d["depth mean"] = other_dp.mean(axis=1)
    qc_d["depth number_1"] = (other_dp > 1).sum(axis=1)
    qc_d["depth number_5"] = (other_dp > 5).sum(axis=1)
    qc_d["depth number_10"] = (other_dp > 10).sum(axis=1)
    qc_d["depth std"] = other_dp.std(axis=1)
    qc_d.columns = [f"{prefix}{x}" for x in qc_d.columns.values]
    return qc_d


def variants_qc(af, dp, groups_d, include_lineage=True):
    """
    :param af: variants-by-cells
    :param dp: variants-by-cells
    :param groups_d: dictionary where group is key and list of cells is value
    :return:
    """
    # Construct vars_qc df for each clone, then merge to make long-df
    all_vars_qc = []
    #all_trans = []
    for group in groups_d:
        print('group', group)
        curr_af = af.loc[:, groups_d[group]]
        curr_dp = dp.loc[:, groups_d[group]]
        qc_d = calc_variants_qc(curr_af, curr_dp,
                         prefix="")
        qc_d["meta"] = "Group"
        qc_d["group"] = group
        qc_d["variant"] = qc_d.index

        other_af = af.loc[:, ~(af.columns.isin(groups_d[group]))]
        other_dp = dp.loc[:, ~(dp.columns.isin(groups_d[group]))]
        other_qc_d = calc_variants_qc(other_af, other_dp,
                                prefix="")
        other_qc_d["variant"] = other_qc_d.index
        other_qc_d["meta"] = "Other"
        other_qc_d["group"] = group
        if include_lineage:
            lin = group.split("_")[0]
            cond = group.split("_")[1]

            qc_d["lineage"] = lin
            other_qc_d["lineage"] = lin
            qc_d["condition"] = cond
            other_qc_d["condition"] = cond
        all_vars_qc.append(qc_d)
        all_vars_qc.append(other_qc_d)

    all_vars_qc_df = pd.concat(all_vars_qc, join='outer', axis=0)
    return all_vars_qc_df, all_vars_qc


def filter_variants(vars_qc, thresh=0.9, thresh_col="perc_0.01"):
    vars_qc = vars_qc.loc[vars_qc[thresh_col]>thresh]
    return vars_qc


def plot_vars(vars_qc, plt_col):
    # 1. Plot variants transitions/transversions for each linege
    sns.countplot(data=vars_qc[vars_qc["meta"] == "Group"], hue="variant type", x="lineage")

    # 2. Plot mean across conditions
    g = sns.catplot(x="condition", y=plt_col, col="lineage", col_wrap=5,
                    data=vars_qc, kind="point", hue="meta",
                    height=2.5, aspect=.8)

    # 3. Plot cell-clone distribution
    sns.countplot(vars_qc["group"])
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


def get_low(df, thresh):
    """ Remove rows that have an average higher than threshold

    :param df:
    :param thresh:
    :return:
    """
    return df.loc[(df.mean(axis=1)<thresh)].index.values


def get_high(df, thresh):
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
