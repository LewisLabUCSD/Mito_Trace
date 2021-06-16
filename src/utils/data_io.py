import pandas as pd
from os.path import join, exists
import os
import logging
import gzip
import os
from scipy.io import mmread
from mplh.fig_utils import helper_save as hs
import pickle

####
## config utils
####
def setup_outdir(outdir, dirs=('figures', 'data','results')):
    for i in dirs:
        if not exists(join(outdir, i)):
            os.mkdir(join(outdir, i))
    return


def save_out(obj, fname, outdir, type, dat_type="pickle", **pd_kwargs):
    """

    :param obj:
    :param outdir:
    :param type:
    :return:
    """
    if type=="figure":
        hs(join(outdir, "figures", fname),f=obj)
    elif type == "data":
        curr_out = join(outdir, "data")
    elif type == "results":
        curr_out = join(outdir, "results")
    else:
        raise ValueError('Please give a valid type (figures, data, results)')

    if dat_type == "pickle":
        pickle.dump(obj, open(join(curr_out, fname.replace('.p', '')+'.p')))
    elif dat_type == "tsv":
        obj.to_csv(join(curr_out, fname.replace('.tsv', '') + '.tsv'), sep='\t', **pd_kwargs)
    elif dat_type == "csv":
        obj.to_csv(join(curr_out, fname.replace('.csv', '') + '.csv'), sep='\t', **pd_kwargs)
    else:
        raise ValueError('Please give a valid dat_type')

    return

####
## df and i/o utils
####
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


def read_csv_multichar(in_f, multicomment=None, encoding='utf-8', **pd_kwargs):
    if multicomment is None:
        return pd.read_csv(in_f, **pd_kwargs)
    else:
        if in_f.endswith(".gz"):
            print('gz')
            f = gzip.open(in_f, 'rt', encoding=encoding)
        else:
            f = open(in_f, "rt")
        curr=0
        for aline in f:
            aline =  str(aline)
            if (len(multicomment) <= len(aline)) and \
                    (multicomment != aline[:len(multicomment)]):
                print(f'skipping {curr} rows')
                break
            else:
                curr += 1
        f.close()
        return pd.read_csv(in_f, skiprows=curr, encoding=encoding,
                           **pd_kwargs)




##############################################
## Vireo spare mmf matrix format utils
##############################################

# Read/write utils
def load_mtx_df(in_f, skip_first=True, give_header=False,
                columns=("Variant", "Cell", "integer")):
    df = pd.read_csv(in_f, comment="%", header=None, sep="\t")
    df.columns = columns
    if skip_first:
        head = df.iloc[0]
        df = df.iloc[1:] # Seems to be summary values
    if give_header and skip_first:
        return df, head
    return df


def wrap_load_mtx_df(indir, oth_f=False, prefix="cellSNP.tag",
                     columns=('Variant', 'Cell', 'integer'), inc_af=False,
                     as_dense=False, var_names=False, vcf_prefix="cellSNP.base"):
    print(prefix)
    curr_ad_f = join(indir, f"{prefix}.AD.mtx")
    curr_dp_f = join(indir, f"{prefix}.DP.mtx")
    if var_names:
        var_meta = read_csv_multichar(join(indir, f"{vcf_prefix}.vcf"),multicomment="##", sep='\t')
        print(var_meta.head())
    if inc_af or as_dense:
        curr_ad = mmread(curr_ad_f).tocsc()
        curr_dp = mmread(curr_dp_f).tocsc()
        AD_df = pd.DataFrame(curr_ad.todense())
        DP_df = pd.DataFrame(curr_dp.todense())
        DP_df.index.name = "Variant"
        DP_df.columns.name = "Cell"
        AD_df.index.name = "Variant"
        AD_df.columns.name = "Cell"
        AD_df = AD_df.sort_index()
        DP_df = DP_df.sort_index()

        if var_names:
            DP_df.index = var_meta.fillna("N").apply(lambda x: str(x["POS"]) + x["REF"]+">"+x["ALT"], axis=1)
            AD_df.index = var_meta.fillna("N").apply(lambda x: str(x["POS"]) + x["REF"]+">"+x["ALT"], axis=1)

        AF_df = AD_df / (DP_df + 0.00001)
        if inc_af:
            return AF_df, DP_df
        else:
            return AD_df, DP_df
    else:
        ad = load_mtx_df(curr_ad_f,
                         columns=columns)
        dp = load_mtx_df(curr_dp_f,
                         columns=columns)
        # if var_names:
        #     ad["Variant"] = ad["Variant"].map(var_meta.index)
        if oth_f:
            oth = load_mtx_df(join(indir, f"{prefix}.OTH.mtx"),
                              columns=columns)
            return ad, dp, oth

    return ad, dp



#def load_with_ids():
#    return


def wrap_write_mtx_df(outdir, ad, dp, oth=None, to_rm=True,
                      prefix="cellSNP.tag", columns=('Variant', 'Cell', 'integer')):

    print(prefix)
    write_mtx_df(ad, outdir, f"{prefix}.AD.mtx", to_rm=to_rm,
                 columns=columns)
    write_mtx_df(dp, outdir, f"{prefix}.DP.mtx", to_rm=to_rm,
                 columns=columns)
    if oth is not None:
        write_mtx_df(oth, outdir, f"{prefix}.OTH.mtx", to_rm=to_rm,
                     columns=columns)
    return


def write_mtx_df(in_mtx, outdir, out_f, to_rm=True, columns=('Variant', 'Cell', 'integer')):
    header = "%%MatrixMarket matrix coordinate integer general\n%\n"
    if to_rm:
        if exists(join(outdir, out_f)):
            os.remove(join(outdir, out_f)) #"cellSNP.tag.DP.mtx"))
    elif exists(join(outdir, out_f)):
        logging.info("File already exists. To rerun, use to_rm")
        return
    with open(join(outdir, out_f), 'a') as file:
        file.write(header)
        full = pd.concat((pd.DataFrame(
            {columns[0]: in_mtx[columns[0]].max(),
             columns[1]: in_mtx[columns[1]].max(),
             columns[2]: in_mtx.shape[0]}, index=["Meta"]),
                             in_mtx.sort_values([columns[0], columns[1]])), sort=False)
        full.to_csv(file, sep="\t", header=False, index=False)
    return


def mgatk_to_vireo(in_af, in_af_meta, in_coverage, outdir, out_name):
    """ Converts mgat output to vireo input.

    Vireo input requires
    :param in_af: tsv file of variant-by-position af
    :param in_af_meta: tsv file for variants (see mgatk). Important columns are position, nucleotide, variant
    :param in_coverage: tsv file like in_af but for coverage
    :param outdir: Output directory to save to
    :param out_name: prefix to the file names

    Output:
    {out_name}.base.vcf.gz: vcf file. Needs #CHROM, POS, REF, ALT columns
    {out_name}.samples.tsv: Cell labels
    {out_name}.tag.AD.mtx vmmf file format
    {out_name}.tag.DP.mtx vmmf file format

    :return:
    """
    #Output: [vcf, AD, DP, cell_labels]
    print('in_coverage', in_coverage)
    print('in_af', in_af)
    af = pd.read_csv(in_af, sep="\t")
    af = af.loc[(~((af==0).all(axis=1)))]
    af = af.loc[:, (~((af == 0).all(axis=0)))]
    af_meta = pd.read_csv(in_af_meta, sep="\t")
    coverage = pd.read_csv(in_coverage, sep="\t")
    af_meta = af_meta.loc[af.index]

    print('af')
    print(af.head())
    print('coverage')
    print(coverage.head())
    coverage = coverage.loc[af.index, af.columns]
    assert(coverage.shape==af.shape)
    assert (coverage.shape[0] == af_meta.shape[0])
    print('coverage')
    print(coverage.head())

    cell_samples = list(af.columns.values)

    # Covert af_meta to vcf-like
    # CHROM  POS     ID_x    REF_x   ALT
    #"position"      "nucleotide"    "variant"

    af_meta["REF"] = af_meta['variant'].apply(lambda x: x.split(">")[0])
    af_meta["ALT"] = af_meta['variant'].apply(
        lambda x: x.split(">")[1])
    af_meta = af_meta.rename({'position':"POS"}, axis=1)
    af_meta["#CHROM"] = "chrM"
    af_meta = af_meta[["#CHROM", "POS", "REF", "ALT",
                       "strand_correlation", "vmr", "n_cells_over_5",
                       "n_cells_over_20"]]

    af_meta = af_meta.sort_values("POS")
    af = af.loc[af_meta.index]
    coverage = coverage.loc[af_meta.index]

    ## Create sparse vmmf files by converting index and columns to 1-based index and then to sparse
    # the order is same for vcf and samples
    import numpy as np
    af.index = np.arange(1, af.shape[0]+1)
    af.columns= np.arange(1, af.shape[1] + 1)
    coverage.index=np.arange(1, coverage.shape[0]+1)
    coverage.columns = np.arange(1, coverage.shape[1] + 1)
    # Convert AF to AD, and then make sparse
    ad = (af*coverage).round().astype(int)
    # Create vmmf files AD and DP
    dp = coverage.rename_axis("Variant").reset_index().melt(id_vars="Variant",
                                      var_name="Cell",
                                      value_name='integer')
    ad = ad.rename_axis("Variant").reset_index().melt(id_vars="Variant",
                                      var_name="Cell",
                                      value_name='integer')

    af_meta.reset_index().to_csv(join(outdir, f"{out_name}.base.vcf"), sep="\t", index=False)
    with open(join(outdir, f"{out_name}.samples.tsv"),'w') as f:
        f.write("\n".join(cell_samples))

    # Dummy oth variable
    oth = pd.DataFrame(columns=ad.columns)
    wrap_write_mtx_df(outdir, ad, dp, oth=oth, to_rm=True, prefix=f"{out_name}.tag",
                      columns=('Variant', 'Cell', 'integer'))
    return
##############################################
## Till here
##############################################
