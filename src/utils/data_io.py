import pandas as pd
from os.path import join, exists, dirname, basename
import logging
import gzip
import os
from scipy.io import mmread
from mplh.fig_utils import helper_save as hs
import pickle
import glob
from icecream import ic
from pandarallel import pandarallel
import tqdm
import numpy as np

####
## config utils
####
def setup_outdir(outdir, dirs=('figures', 'data','results')):
    if not exists(outdir):
        print("creating outdir", outdir)
        os.mkdir(outdir)
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


def import_from_dir(indir, prefix="", suffix="", f_list=None):
    if prefix is None:
        prefix = ""
    if f_list is not None:
        files = []
        for i in f_list:
            curr_f = glob.glob(os.path.join(indir,
                                            prefix+i+suffix))
            if len(curr_f) > 1:
                print(f"Two files show up for {i}. Please give more specific strings to specify the files")
                return []
            elif len(curr_f) == 0:
                print(f"{i} not found in directory {indir}. Continuing "
                      f"without adding it")
            else:
                files.append(curr_f[0])
    else:
        files = glob.glob(os.path.join(indir, prefix+"*"+suffix))
    return files


####
## df and i/o utils
####
def combine_dfs(dfs, use_key=True, key_label="Key", na_val=0,
                update_col_names=False, col_inds=None, axis=1,
                return_meta=False):
    """

    :param dfs: dict of dfs, where keys are names to be used for a new column,
                and/or for add that as a suffix to the column names.
    """
    meta = {}
    for key in dfs:
        if use_key:
            dfs[key][key_label] = key
            meta[key] = dfs[key][key_label]
        if update_col_names:
            dfs[key].columns = dfs[key].columns.astype(str) + "_" + str(key)
        if col_inds is not None:
            dfs[key] = dfs[key].loc[:, dfs[key].columns.isin(col_inds)]

    if return_meta:
        all = pd.concat(dfs.values(), axis=axis, sort=False).fillna(na_val)
        if use_key:
            all = all.drop(key_label, axis=1)
        return all, dfs

    return pd.concat(dfs.values(), axis=axis, sort=False).fillna(
        na_val)  # return  pd.concat(dfs.values()).fillna(fillna)


def read_csv_multichar(in_f, multicomment=None, encoding='utf-8', verbose=False, **pd_kwargs):
    if verbose:
        ic.enable()
    else:
        ic.disable()
    if multicomment is None:
        return pd.read_csv(in_f, **pd_kwargs)
    else:
        if in_f.endswith(".gz"):
            ic('gz')
            f = gzip.open(in_f, 'rt', encoding=encoding)
        else:
            f = open(in_f, "rt")
        curr=0
        for aline in f:
            aline =  str(aline)
            if (len(multicomment) <= len(aline)) and \
                    (multicomment != aline[:len(multicomment)]):
                ic(f'skipping {curr} rows')
                break
            else:
                curr += 1
        f.close()
        return pd.read_csv(in_f, skiprows=curr, encoding=encoding,
                           **pd_kwargs)




##############################################
## single-cell pileup data for each nucleotide. 'scPileup' format.
##############################################
def load_sc_pileup(sc_f):
    """ Loads and sets the columns of the coverage pileup.
    :param sc_f:
    :return:
    """
    sc_coverage = pd.read_csv(glob.glob(sc_f)[0], header=None)
    print(sc_coverage.head())
    if len(sc_coverage.columns)>3:
        sc_coverage.columns = ["Position", "CB", "Coverage", "Count Fw", "Count Rev"]
        sc_coverage["Count Fw"] = sc_coverage["Count Fw"].astype(np.int64)
        sc_coverage["Count Rev"] = sc_coverage["Count Rev"].astype(np.int64)
    else:
        sc_coverage.columns = ["Position", "Cell", "Coverage"]
    sc_coverage["Cell"] = sc_coverage["Cell"].apply(lambda x: x.replace(".bam",""))
    sc_coverage["Coverage"] = sc_coverage["Coverage"].astype(np.int64)
    logging.info('sc_coverage head - {}'.format(sc_coverage.head().to_string()))
    return sc_coverage


def load_and_filter_nt_pileups(concat_dir, cell_inds, nt_pos,
                               out_d=None, name=None, incl_cov=False,
                               sum_cov=False, input_suffix=".strands.txt.gz",
                               compression=False):

    """ Filters the NT pileup and bq matrices and saves as pileups.

    :param concat_dir:
    :param cell_inds:
    :param nt_pos:
    :param out_d:
    :param name:
    :return:
    """
    print("Loading pileups")
    nt_pileup = {}
    for n in ["A", "C", "G", "T"]:
        curr_f = glob.glob(os.path.join(concat_dir, f"*.{n}{input_suffix}"))[0]
        df = pd.read_csv(curr_f, header=None)
        # Only take the Forward reads since it is the one to compare against the reference
        if len(df.columns) > 4:
            print('num cols', len(df.columns))
            df.columns = ["Position", "Cell", "Fw Coverage", "Fw BQ",
                          "Rev Coverage", "Rev BQ"]
            df = df.fillna(0)
            df["Fw Coverage"] = df["Fw Coverage"].copy().astype(np.int64)
            df["Rev Coverage"] = df["Rev Coverage"].copy().astype(np.int64)
            curr_cov = df["Fw Coverage"] + df["Rev Coverage"]
            df["BQ"] = df["Fw BQ"]*(df["Fw Coverage"]/curr_cov) + df["Rev BQ"] * (df["Rev Coverage"]/curr_cov)
            if sum_cov:
                df["Coverage"] = curr_cov
        else:
            df.columns = ["Position", "Cell", "Coverage", "BQ"]
        df["Cell"] = df["Cell"].apply(lambda x: x.replace(".bam", ""))
        df = df.copy()[df["Cell"].isin(cell_inds)]
        df = df.copy()[df["Position"].isin(nt_pos)]
        if not (out_d is None or name is None):
            curr_out = join(out_d, f"{name}.{n}.txt")
            #Drop any 0s
            if "Coverage" in df.columns:
                df = df.drop("Coverage", axis=1)
            if "BQ" in df.columns:
                df = df.drop("BQ", axis=1)
            df = df.sort_values(["Cell", "Position"])
            print('df no 0 shape', df[(df["Fw Coverage"] !=0 ) & (df["Rev Coverage"]!=0)].shape)
            df[(df["Fw Coverage"] !=0 ) & (df["Rev Coverage"]!=0)].to_csv(curr_out, header=None, index=False)
        nt_pileup[n] = df

    if incl_cov:
        df = load_sc_pileup(
            os.path.join(concat_dir, f"*.coverage{input_suffix}"))
        df = df.copy()[df["Cell"].isin(cell_inds)]
        df = df.copy()[df["Position"].isin(nt_pos)]
        df = df.fillna(0)

        if not (out_d is None or name is None):
            curr_out = join(out_d, f"{name}.coverage.txt")
            df = df.sort_values(["Cell", "Position"])
            df[(df["Coverage"] != 0)].to_csv(curr_out, header=None, index=False)
            df.to_csv(curr_out, header=None, index=False)
        nt_pileup["coverage"] = df
    return nt_pileup


##############################################
## Vireo spare mmf matrix format utils
##############################################
# Read/write utils
def load_mtx_df(in_f, skip_first=True, give_header=False,
                columns=("Variant", "Cell", "integer"), sep="\t"):
    df = pd.read_csv(in_f, comment="%", header=None, sep=sep)
    df.columns = columns
    if skip_first:
        head = df.iloc[0]
        df = df.iloc[1:] # Number of elements listed here
    if give_header and skip_first:
        return df, head
    return df


def wrap_load_mtx_df(indir, oth_f=False, prefix="cellSNP.tag",
                     columns=('Variant', 'Cell', 'integer'), inc_af=False,
                     as_dense=False, var_names=False, vcf_prefix="cellSNP.base",
                     cell_names=False, cell_prefix="cellSNP.samples.tsv",
                     verbose=True):
    if not verbose:
        ic.disable()
    else:
        ic.enable()
    ic(prefix)
    curr_ad_f = join(indir, f"{prefix}.AD.mtx")
    curr_dp_f = join(indir, f"{prefix}.DP.mtx")
    if var_names:
        var_meta = read_csv_multichar(join(indir, f"{vcf_prefix}.vcf"),multicomment="##", sep='\t')
        ic(var_meta.head())
    if cell_names:
        cells = pd.read_csv(join(indir, cell_prefix), header=None)[0].values
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

        if cell_names:
            AD_df.columns = cells
            DP_df.columns = cells

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
def write_mtx_df(in_mtx, outdir, out_f, to_rm=True,
                 columns=('Variant', 'Cell', 'integer'), sep="\t"):
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
        full.to_csv(file, sep=sep, header=False, index=False)
    return


def wrap_write_mtx_df(outdir, ad, dp, oth=None, to_rm=True,
                      prefix="cellSNP.tag",
                      columns=('Variant', 'Cell', 'integer')):

    print(prefix)
    write_mtx_df(ad, outdir, f"{prefix}.AD.mtx", to_rm=to_rm,
                 columns=columns)
    write_mtx_df(dp, outdir, f"{prefix}.DP.mtx", to_rm=to_rm,
                 columns=columns)
    if oth is not None:
        write_mtx_df(oth, outdir, f"{prefix}.OTH.mtx", to_rm=to_rm,
                     columns=columns)
    return




###### MGATK output converted to the Vireo specified input ######
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
    # Remove cells and positions that are all 0s
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
    af_meta = af_meta.sort_values("POS")

    af = af.loc[af_meta.index]
    coverage = coverage.loc[af_meta.index]

    ## Create sparse vmmf files by converting index and columns to 1-based index and then to sparse
    # the order is same for vcf and samples

    af.index = np.arange(1, af.shape[0]+1)
    af.columns= np.arange(1, af.shape[1] + 1)
    coverage.index = np.arange(1, coverage.shape[0]+1)
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

    # Save files
    af_meta = af_meta.reset_index()
    af_meta = af_meta[["#CHROM", "POS", "REF", "ALT",
                       "strand_correlation", "vmr", "n_cells_over_5",
                       "n_cells_over_20", "index"]]
    af_meta.to_csv(join(outdir, f"{out_name}.base.vcf"), sep="\t", index=False)

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


##############################################
## Filtered output of cellranger -> fragments.tsv.gz
##############################################
def sparse_to_cellranger_fragments(mtx_f, peaks_f, cells_f, out_f, cell_col="Group_Barcode", n_cpus=2):
    """ Filtered output of cellranger -> fragments.tsv.gz

    :param mtx_f:
    :param peaks_f:
    :param cells_f:
    :param out_f:
    :param cell_col:
    :param n_cpus:
    :return:
    """
    ic(mtx_f, peaks_f, cells_f, out_f, cell_col)
    frags = load_mtx_df(mtx_f, sep=' ')
    #frags = frags.pivot(index='Variant', columns='Cell', values='integer')
    #frags = frags.sort_index().sort_index(axis=1)
    peaks = pd.read_csv(peaks_f, sep='\t', header=None, comment="#")
    peaks.index = peaks.index + 1
    peaks = peaks.rename({0:"Chr", 1:"Start", 2:"End"}, axis=1)
    cells = pd.read_csv(cells_f, sep='\t')[cell_col]
    cells.index = cells.index+1
    frags["CB"] = frags["Cell"].map(cells)
    frags = pd.merge(frags, peaks, left_on="Variant", right_index=True)
    frags = frags.sort_values(["Chr", "Start", "End"])
    # Add bed file information (This is the information that is lost btwn using fragments and the bed file.
    frags[['Chr','Start','End', "CB", "integer"]].to_csv(out_f, sep='\t', index=False, header= False, compression='gzip')

    pandarallel.initialize(nb_workers=n_cpus)
    # Create position indices with format "chr:start-end"
    frags["ID"] = frags.parallel_apply(lambda x: f"{x['Chr']}:{x['Start']}-{x['End']}", axis=1)\

    # Create sparse dgcMatrix, with index as peak, column CB, value is integer
    frags[["ID", "CB", "integer"]].to_csv(f"{out_f}.counts.mtx",
                                          sep="\t", header=None,
                                          index=False)
    return frags


def add_id_to_sparse(sparse, cell=None, pos=None,
                     cell_col=None, pos_col=None,
                     cell_name="Cell", pos_name="Pos", drop=False):
    """ Adds id's to the indexed matrices using the index maps for cells
    and the values

    :param sparse: Sparse df, with  at least columns cell_col and pos_col.
    :param cell: Map of index:cell_barcode where index is 1-based to the matrix file
    :param pos: Map of index:data_barcode
    :param cell_col: Column for the cell
    :param pos_col: Column for data
    :param cell_name: New name to use for cell (default: "Cell")
    :param pos_name: New name to use for position (default: "Pos")
    :param drop: To drop the iniial columns
    :return:
    """
    if cell is not None:
        sparse[cell_name] = sparse[cell_col].map(cell)
    if pos is not None:
        sparse[pos_name] = sparse[pos_col].map(pos)

    if cell is not None and (drop and cell_col!=cell_name):
        sparse = sparse.drop([cell_col], axis=1)
    if pos is not None and (drop and pos_col != pos_name):
        sparse = sparse.drop([pos_col], axis=1)
    return sparse


def get_filt_indir(allConfig, filt_prefix, covDir_dict, use_cov_f=False):
    allFilt = {}
    for c in allConfig:
        #print(c)
        for s in allConfig[c]["multiplex"]["samples"]:
            ###########################
            # Get filtered cell indices
            ###########################
            f_in_data = join(allConfig[c]["results"], "data", s, "MT",
                             "cellr_True", f"{s}_200", filt_prefix,
                             "af_by_cell.DP.tsv")
            f_in = join(allConfig[c]["results"], s, filt_prefix,
                        "af_by_cell.DP.tsv")
            if "old_results" in allConfig[c]:
                f_in_old = join(allConfig[c]["old_results"], s, "MT",
                                "cellr_True", f"{s}_200", filt_prefix,
                                "af_by_cell.DP.tsv")
            else:
                f_in_old = ""
            if c in covDir_dict:
                f_in_covDir_dict = join(covDir_dict[c], s, filt_prefix,
                                        "af_by_cell.DP.tsv")
            else:
                f_in_covDir_dict = ""
            does_exist = False
            for f in [f_in, f_in_data, f_in_covDir_dict, f_in_old]:
                if exists(f):
                    #print(f"Adding {f}")
                    if use_cov_f:
                        allFilt[(c, s)] = f.replace("af_by_cell.DP.tsv",
                                                    f"{s}.coverage.txt")
                    else:
                        allFilt[(c, s)] = f
                    does_exist = True
            if not does_exist:
                print(f"File for {c}, {s} not here")
                continue
    return allFilt


def get_cov_indir(allConfig, covDir_dict, input_suffix=".strands.txt.gz"):
    allCov = {}
    for c in allConfig:
        for s in allConfig[c]["multiplex"]["samples"]:
            f_in_data = join(allConfig[c]["results"], "data", s,
                                 "MT","cellr_True", f"{s}_200",
                                 f"{s}.coverage{input_suffix}")
            f_in = join(allConfig[c]["results"], s, "MT","cellr_True", f"{s}_200",
                                     f"{s}.coverage{input_suffix}")
            if "old_results" in allConfig[c]:
                f_in_old = join(allConfig[c]["old_results"], s,
                                         "MT","cellr_True", f"{s}_200",
                                         f"{s}.coverage{input_suffix}")
            else:
                f_in_old = ""
            if c in covDir_dict:
                f_in_covDir_dict = join(covDir_dict[c], s,
                                         "MT","cellr_True", f"{s}_200",
                                         f"{s}.coverage{input_suffix}")
            else:
                f_in_covDir_dict = ""
            does_exist=False
            for f in [f_in, f_in_data, f_in_covDir_dict, f_in_old]:
                if exists(f):
                    print(f"Using {f}")
                    does_exist=True
                    allCov[(c,s)] = f
            if not does_exist:
                print(f"File for {c}, {s} not here")
                continue
    return allCov


def preproc_cellr_indir(allConfig, c, s):
    samples = pd.read_table(allConfig[c]["samples"], sep=',',
                            index_col=0).reset_index().set_index(
        "sample_name", drop=False)
    samples = samples.dropna(axis=1)
    if s in samples.index:
        curr_in = dirname(samples.loc[s, "barcode_f"])
    else:
        curr_in = dirname(
            samples.set_index("sample").loc[s, "barcode_f"])
    print(c, s)
    curr_frags = pd.read_csv(join(dirname(curr_in), "fragments.tsv.gz"),
                             sep="\t", header=None)
    curr_frags.columns = ["Chr", "Start", "End", "Cell", "Count"]
    curr_frags["Cell"] = curr_frags["Cell"] + "_" + s + "_" + c
    peaks = pd.read_csv(join(curr_in, "peaks.bed"), sep="\t",
                        header=None)
    peaks = peaks.rename({0: "Chr", 1: "Start", 2: "End"}, axis=1)

    cells = pd.read_csv(join(curr_in, "barcodes.tsv"), header=None)[0]
    cells = cells + "_" + s + "_" + c
    cells.index += 1  # 1-based map

    frags = {"Cells": cells, "Peaks": peaks, "Frags": curr_frags}
    return frags



#def read_
# if __name__ == "__main__":
#     sparse_to_cellranger_fragments("/data2/mito_lineage/data/external/granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.mtx",
#                                    "/data2/mito_lineage/data/external/granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.peaks.bed",
#                                    "/data2/mito_lineage/data/external/granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.cell_barcodes.txt",
#                                    out_f="/data2/mito_lineage/data/processed/external/granja_cd34/granja_cd34.fragments.tsv")
