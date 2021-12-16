import src.utils.variant_utils as vu
import pandas as pd
import seaborn as sns
from mplh.fig_utils import helper_save as hs
from os.path import join


def calc_fold(af, curr_cells_meta, samples, vars_to_plot):
    variant_df = vu.type_of_variants(af.index)

    variants_box = vu.variants_dense(af, vars_to_plot, samples_d=curr_cells_meta['condition'],
                                     donors_d=curr_cells_meta['donor'],
                                     variant_d=variant_df['variant type'])


    fold_vars = {}
    mean_vars = {}
    for ind, val in variants_box.groupby("Variant"):
        cond_means = (val.groupby("condition").mean())[["sqrtAF"]]
        mean_vars[ind] = cond_means.values.flatten()

        cond_fold = cond_means.loc[samples[1]] / cond_means.loc[samples[0]]
        fold_vars[ind] = cond_fold

    fold_vars = pd.DataFrame(fold_vars)

    len(variants_box["Variant"].unique())

    mean_vars_df = pd.DataFrame(mean_vars, index=samples)
    mean_vars_df = mean_vars_df.reset_index().rename({"index": "condition"}, axis=1).melt(id_vars="condition")

    return pd.DataFrame(fold_vars), mean_vars_df


def cells_thresh_reads(var_r_rt, t, rt):
    return var_r_rt[(var_r_rt["af"] > t) & (var_r_rt["reads"] > rt)].index

def cells_thresh(var_series, thresh):
    return var_series[var_series > thresh].index


def cells_median(var_series, cells):
    return var_series.loc[cells].median()


def cells_mean(var_series, cells):
    return var_series.loc[cells].mean()


def calc_variant_descriptives(af, vars_df, prefix, thresholds, dp, read_thresh):
    vars_cells = {}
    for t in thresholds:
        vars_cells[t] = af.apply(cells_thresh, args=(t,))
        # for rt in read_thresh:
        #     curr_cells = af.apply(cells_thresh_reads, args=(t, rt))
        #     vars_cells[(t, rt)] = af.loc[curr_cells].apply(cells_thresh_reads, args=(t, rt))
    vars_df[f"{prefix}median"] = af.mean(axis=0)
    vars_df[f"{prefix}mean"] = af.mean(axis=0)

    vars_df[f"{prefix}coverage_median"] = dp.mean(axis=0)
    vars_df[f"{prefix}coverage_mean"] = dp.mean(axis=0)

    for t in thresholds:
        vars_df[f"{prefix}nCells_{t}"] = len(vars_cells[t])
        vars_df[f"{prefix}median_{t}"] = af.apply(cells_median, args=(vars_cells[t],))
        # for rt in read_thresh:
        #     vars_df[f"{prefix}nCells_{t}_{rt}"] = len(vars_cells[(t, rt)])
            #vars_df[f"{prefix}median_{t}_{rt}"] = dp.apply(cells_median, args=(vars_cells[t],))
    return vars_df


def vars_as_clones(af, dp, cells_meta, samples, thresholds, read_thresh, groups=["condition"]):
    """ Creates cell clones based on single MT variants

    af: cell-by-variant af, where the cell has a unique id, and the variants are in format '{REF}{POS}>{ALT}'
    cells_meta: Ways to group the cells, with a 'condition' column being used
    """

    # Steps:
    # 1. create variant df with the proper columns
    # 2. get cell ID for having the variant, > thresh for each. Add in NA if not passing read filter.
    #    get cell IDs for not having the variant,
    # 3. get #cells in condition > thresh for each.
    # 4. get #cells > thresh for each
    # 5. calculate median, mean for the cells that have the variant.
    # 6. calculate median, mean overall
    # 7. Plot distribution of variant AF , only for ones with > 10 cells.
    # 8. Correlate # cells with variant and fold change
    # 9. Calculate nearest shared variant for each by calculating #overlap for each

    vars_df = vu.type_of_variants(af.columns, style=">", to_preproc=True)
    vars_df = calc_variant_descriptives(af, vars_df, prefix="", thresholds=thresholds, dp=dp, read_thresh=read_thresh)

    if cells_meta is not None:
        for c, df in cells_meta.groupby(groups):
            vars_df = calc_variant_descriptives(af.loc[df.index], vars_df,
                                                prefix=f"condition {c}_",
                                                thresholds=thresholds, read_thresh=read_thresh)

        # fold change cell # and AF
        for t in thresholds:
            vars_df[f"FoldChange_nCells_{t}"] = (vars_df[f"condition {samples[1]}_nCells_{t}"])/(vars_df[f"condition {samples[0]}_nCells_{t}"])
            # plot FC vs number of cells total
            vars_df[f"FoldChange_medianAF_{t}"] = (vars_df[f"condition {samples[1]}_median"]) / (vars_df[f"condition {samples[0]}_median_{t}"])

        vars_df[f"FoldChange_medianAFall"] = (vars_df[f"condition {samples[1]}_median"]) / (vars_df[f"condition {samples[0]}_median"])
    return vars_df


def plots_vars(vars_df, thresholds, name, outdir):
    for t in thresholds:
        sns.regplot(vars_df, x=f"FoldChange_nCells_{t}", y=f"FoldChange_medianAF_{t}")
        hs(join(outdir, f"{name}_t"), to_pdf=False)
        # plot barplot across conditions
        # for the ordered sizes, see if transition/transversion, see if other cells have coverage at that position,
        # compare to other donors, if we have strand concordance plot that too
    return

def save_clones(vars_dict):
    return


def main(af_f, dp_f, cells_meta_f, samples, name, outdir, clone_thresh=10):
    af = pd.read_csv(af_f, sep="\t")
    dp = pd.read_csv(dp_f, sep="\t")

    cells_meta = pd.read_csv(cells_meta_f)
    samples = samples.split(",")
    thresholds = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9]
    read_thresh = [2, 10, 100]
    vars_df = vars_as_clones(af, cells_meta, samples, thresholds, read_thresh, groups=["condition"])
    plots_vars(vars_df, thresholds, name, outdir)
    # save ones with cell_filter count

    return
