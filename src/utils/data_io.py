import pandas as pd
from os.path import join

def combine_dfs(dfs, use_key=True, key_label="Key", naval=0,
                update_col_names=False, col_inds=None, axis=1):
    """

    :param dfs: dict of dfs, where keys are names to be used for a new column,
                and/or for add that as a suffix to the column names.
    """
    for key in dfs:
        if use_key:
            dfs[key][key_label] = key
        if update_col_names:
            dfs[key].columns = dfs[key].columns + "_" + key
        if col_inds is not None:
            dfs[key] = dfs[key].loc[:, dfs[key].columns.isin(col_inds)]
    return pd.concat(dfs.values(), axis=axis, sort=False).fillna(
        naval)  # return  pd.concat(dfs.values()).fillna(fillna)


def dense_to_sparse():
    # Output: [row, col, val] Keys
    return


def sparse_to_dense():
    return


def mgatk_to_vireo(in_af, in_af_meta, in_coverage, outdir, out_name):
    #Output: [vcf, AD, DP, cell_labels]
    af = pd.read_csv(in_af, sep="\t")
    af_meta = pd.read_csv(in_af_meta, sep="\t")
    coverage = pd.read_csv(in_coverage, sep="\t", header=None)

    coverage = coverage[coverage[1].isin(af.columns)] #af_meta["nucleotide"]
    coverage = coverage[coverage[0].isin(af_meta['"Position"'])]

    # coverage can now be saved.
    coverage.to_csv(join(outdir, f"{out_name}.DP.txt" ))

    # Convert AF to AD, and then make sparse
    af = af*(coverage.pivot(index=0, columns=1, values=2))
    print("af")
    print(af.head())
    af.to_csv(join(outdir, f"{out_name}.AD.txt"))

    # Covert af_meta to vcf-like
    # CHROM  POS     ID_x    REF_x   ALT
    #"position"      "nucleotide"    "variant"
    af_meta["REF"] = af_meta['"variant'].apply(lambda x: x.split(">")[0])
    af_meta["ALT"] = af_meta['"variant'].apply(
        lambda x: x.split(">")[1])
    af_meta.rename({"position":"POS"}, axis=1)
    af_meta["#CHROM"] = "chrM"
    af_meta = af_meta[["#CHROM", "POS", "REF", "ALT",
                       "strand_correlation", "vmr", "n_cells_over_5",
                       "n_cells_over_20"]]
    af_meta.to_csv(join(outdir, f"cellSNP.base.vcf"))
    return
