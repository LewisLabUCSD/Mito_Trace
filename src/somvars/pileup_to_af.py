import pandas as pd
from pandarallel import pandarallel
#pandarallel.initialize(16)


def is_var(entry, vcf_df):
    if entry["ID"] in vcf_df.index:
        return True
    return False


def get_proper_variants(df, vcf_df):
    df["ID"] = df.parallel_apply(
        lambda x: f"{x['posID']}_{x['nt']}", axis=1)
    df["isVar"] = df.parallel_apply(is_var, args=(vcf_df,), axis=1)
    return df


def get_counts(df):
    df["coverage"] = df.groupby(["chr", "pos", "cell"])["count"].transform(lambda x: sum(x))
    df["AF"] = df["count"] / df["coverage"]
    return df


def filt_pos(df, vcf_df):
    pileup_vars_df = df.loc[
        df["posID"].isin(vcf_df["posID"])].copy()

    vcf_pos_df = vcf_df.set_index("posID")
    vcf_pos_df = vcf_pos_df.loc[~(
        vcf_pos_df.index.duplicated())]  # rm duplicated position variants just for position
    vcf_pos_df.head()

    pileup_vars_df["ref"] = pileup_vars_df.apply(
        lambda x: vcf_pos_df.loc[x["posID"], "REF"], axis=1)

    pileup_vars_df = df.loc[
        df["posID"].isin(vcf_df["posID"])].copy()
    return pileup_vars_df


def create_pileup(pileup_df):
    #pileup_df = pd.read_csv(chunk, sep="\t", index_col=0, compression='gzip', dtype=object)

    ### When concatenating the column names were added too so need to drop those entries and convert the dtypes
    print(pileup_df.shape)
    pileup_df = pileup_df.loc[~((pileup_df["chr"] == "chr") & (pileup_df["nt"] == "nt"))].copy()
    # Add posID and BQ and change to proper dtypes
    pileup_df["posID"] = pileup_df["chr"] + "_" + pileup_df["pos"]
    pileup_df = pileup_df.astype({"pos":int, "BQ": float, "count": int})
    pileup_df["BQ"] = pileup_df["BQ"].astype(int)
    return pileup_df


def process_chunk(chunk, vcf_df):
    """ Gets the pileup_vars_df for each chunk"""
    pileup_df = create_pileup(chunk)
    pileup_vars_df = filt_pos(pileup_df, vcf_df)
    pileup_vars_df = get_counts(pileup_vars_df)
    vcf_df["var"] = vcf_df["posID"] + "_" + vcf_df["ALT"]
    vcf_df = vcf_df.set_index("var")
    pileup_vars_df = get_proper_variants(pileup_vars_df, vcf_df)
    return pileup_vars_df



