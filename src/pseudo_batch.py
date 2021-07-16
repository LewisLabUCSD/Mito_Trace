""" Merges vireo inputs together across different conditions.

## Input:
outdir, type=click.Path()): The directory to save output.
indirs, type=click.Path(exists=True): Input directories.
--cell_subsamp:
--is_prop, default=False, type=click.BOOL: If True, subsampling preserves the proportion of each conditions cell counts.
--num_cells, default=1000: Number of cells total to sample. If above the total number of cells, will use all cells.
--samples, default="", Comma-separated string of the condition names. This will be added to the suffix of the cell labels.


### Required files in input directories:
- cellSNP.base.vcf.gz: tsv file with required columns ["#CHROM", "POS", "REF", "ALT"]
- cellSNP.samples.tsv: Line separated file of cell labels.

Vireo Matrix Market Format (vmmf)
Matrix market format files with the first two lines are comments (contain '%').
Sparse matrix (1-based index), where first column is position index, 2nd column is cell index, and 3rd is the count.
The third line also has 3 columns, and the values are a) # of unique indices in 1st column,
                                                      b) # indices in 2nd column,
                                                      c) total number of observations (i.e. number of rows)
Cell index refers to the cellSNP.samples.tsv file.
Position index refers to cellSNP.base.vcf.gz file

- cellSNP.tag.AD.mtx: Alternative allele depth in matrix market format (see below)
- cellSNP.tag.DP.mtx: Total Depth in martix market format.
- (optional) cellSNP.tag.OTH.mtx: Other (non Alt or Ref allele) depth.

## Output files:
1. variant_indices_{sampleIndex}.tsv: tsv file for each initial sample (index corresponds to the order of the input)
                                  2 columns with no header, where first column is 1-based index of initial input vcf for that condition,
                                  and the second column is the 1-based index for the output merged vcf file.
2. cell_indices_{sampleIndex}.txt: csv file column headers a) old index, b) new index. Also 1-based index

3. cell_labels.txt: csv file of cell labels with sample ID added as suffix.
        header: ID,raw ID,new index
        row example: AAACGAATCTTACTCA-1_Control,AAACGAATCTTACTCA-1,1
4. cellSNP.base.vcf: vcf file as described above.
5. cellSNP.tag.AD.mtx: same as above, vmmf files.
6. cellSNP.tag.DP.mtx


## Temporary classes/variables used:
vmmf_df: the sparse matrix, columns are Variant, Cell, integer
"""
from os.path import join
import os
from numpy import random
import numpy as np
import pandas as pd
import click
import sys
import logging
from src.utils.data_io import wrap_write_mtx_df, load_mtx_df, wrap_load_mtx_df, read_csv_multichar


def run_batch(indirs, outdir, num_reads_total):
    num_reads_each = len(indirs)/num_reads_total
    return


########################################
def add_suffix_to_labels(in_files, cells_kept=None, out_file=None, sample_names=None):
    """
    :param in_files: list of labels files
    :param out_file: output label file
    :return: df of the filtered cell list.

    columns are 'ID', 'raw ID', and 'new index'. The first contains the
    suffix with the old id, the new index contains the mapping to the
    outputted subsampled cells, which is 1-based and the raw ID is the initial cell IDs without the suffix.
    """
    print('sample_names', sample_names)
    all = []
    for ind, f in enumerate(in_files):
        curr = pd.read_csv(f, header=None)
        curr["raw ID"] = curr[0]
        if (sample_names is None) or (sample_names==[]):
            curr[0] = curr[0] + "_" + os.path.basename(f)
        else:
            curr[0] = curr[0] + "_" + sample_names[ind]
        if cells_kept is not None: # Indices to keep
            inds_to_keep = np.array(list(cells_kept[os.path.dirname(f)].keys()))-1
            curr = curr.iloc[inds_to_keep]
            curr['new index'] = list(cells_kept[os.path.dirname(f)].values())
            curr["condition"] = sample_names[ind]
        curr = curr.rename({0: "ID"}, axis=1)
        all.append(curr)
    df = pd.concat(all, ignore_index=True)

    if out_file is not None:
        df.to_csv(out_file, index=False, sep='\t')
    return df
########################################


def merge_vcf_ids(indirs, outdir=None,prefix_vcf_in="cellSNP.base",
                  prefix_vcf_out="cellSNP.base"):
    """ Merges the variant files called from cellSNP

    :param indirs:
    :return:
    """

    curr_vcf = join(indirs[0], f"{prefix_vcf_in}.vcf")
    if not os.path.exists(curr_vcf):
        if not os.path.exists(curr_vcf + ".gz"):
            raise ValueError(f"VCF file {curr_vcf} not here!")
        else:
            curr_vcf = curr_vcf + ".gz"

    variants = read_csv_multichar(curr_vcf, multicomment="##", sep='\t') # pd.read_csv(curr_vcf, skiprows=1, sep="\t")
    #variants = pd.read_csv(curr_vcf, skiprows=1, sep="\t")
    variants["old " + indirs[0]] = variants.index.values.astype(int) + 1

    old_variants = {}
    old_variants[indirs[0]] = variants.copy()
    old_variants[indirs[0]].loc[:, "old"] = 0
    old_variants[indirs[0]].loc[:, "old"] = old_variants[
                                                indirs[0]].index + 1

    for ind, val in enumerate(indirs[1:]):
        curr_vcf = join(val, f"{prefix_vcf_in}.vcf")
        if not os.path.exists(curr_vcf):
            if not os.path.exists(curr_vcf + ".gz"):
                raise ValueError(f"VCF file {curr_vcf} not here!")
            else:
                curr_vcf = curr_vcf + ".gz"
        curr = read_csv_multichar(curr_vcf, multicomment="##", sep='\t')
        # curr["old index"] = curr.index.values+1
        # curr[val] = curr.index.values + 1
        curr["old " + val] = curr.index.values.astype(int) + 1
        variants = pd.merge(variants, curr, on=["index", "#CHROM", "POS", "REF", "ALT"],
                            how="outer")

        old_variants[val] = curr.copy()
        old_variants[val]["old"] = 0
        old_variants[val]["old"] = curr.index.values + 1

    print('variants')
    print(variants.head())
    print(variants.tail())

    # Loop again and map the coordinates
    vars_coords = dict()
    variants["new ID"] = variants.index.values + 1
    full_vars = variants.copy()

    # what the new index should be.
    for val in indirs:
        curr_vcf = join(val, f"{prefix_vcf_in}.vcf")
        if not os.path.exists(curr_vcf):
            if not os.path.exists(curr_vcf + ".gz"):
                raise ValueError(f"VCF file {curr_vcf} not here!")
            else:
                curr_vcf = curr_vcf + ".gz"
        curr_vars = read_csv_multichar(curr_vcf, multicomment="##", sep='\t') # pd.read_csv(curr_vcf, skiprows=1, sep="\t")
        curr_vars["old index"] = curr_vars.index + 1
        curr_vars = pd.merge(curr_vars, full_vars, how="inner",
                             on=["#CHROM", "POS", "ALT"])
        curr_vars.index = curr_vars["old index"]
        vars_coords[val] = curr_vars["new ID"]  # Pandas series

    if outdir is not None:
        full_vars.to_csv(join(outdir, f"{prefix_vcf_out}.vcf"), sep='\t',
                        index=False)
        for ind, val in enumerate(indirs):
            out_f = join(outdir, f"variant_indices_{ind}.tsv")
            vars_coords[val].index.name = val
            vars_coords[val].to_csv(out_f,sep="\t")  # Use index + header
    return vars_coords, full_vars


def subsample_sparse_matrices(outdir, indirs, cell_subsample=0.1,
                              num_cells_total=1000,
                              is_proportional=False, sample_names=None,
                              prefix_tags_in="cellSNP.tag",
                              prefix_vcf_in="cellSNP.base",
                              prefix_tags_out="cellSNP.tag",
                              prefix_vcf_out="cellSNP.base",
                              prefix_samples_in="cellSNP"):
    """ Subsamples the cellSNP output matrices and combines into a new folder, with metadata saying which cell comes from which sample.

    :param indirs:
    :param outdir:
    :param cell_subsample:
    :param num_cells_total:
    :param is_proportional:
    :return:
    """
    print(outdir, indirs, cell_subsample, num_cells_total,
          is_proportional)
    print(type(is_proportional))
    if outdir in indirs:
        raise ValueError(
            "Outdir is one of the indirs. This cant be. Please change one.")

    # 1. Merge VCF files (union of them) and save new indices
    vars_coords, variants = merge_vcf_ids(indirs, outdir=outdir,
                                          prefix_vcf_in=prefix_vcf_in,
                                          prefix_vcf_out=prefix_vcf_out)
    # 2. Count the number of cells in each input to determine proportions
    cell_count = {}
    for i in indirs:
        ad, dp, oth = wrap_load_mtx_df(i, oth_f=True, prefix=prefix_tags_in,
                     columns=('Variant', 'Cell', 'integer'))
        if len(ad) == 0:
            raise ValueError("The allele matrix is empty. It could have all been filtered out.")
        # cell_count[i] = len(set(dp["Cell"].values))
        cell_count[i] = len(list(
            set(dp["Cell"].values).union(set(ad["Cell"].values)).union(
                set(oth["Cell"].values))))
    cell_proportion = np.array(list(cell_count.values())) / sum(
        cell_count.values())
    print('cell proportion', cell_proportion)

    # ad["Variant"] = vars_coords[ad["Variant"]]
    # dp["Variant"] = vars_coords[dp["Variant"]]
    # oth["Variant"] = vars_coords[oth["Variant"]]

    # 3. For each input, subsample based on either proportion or cell number, and remap the indices.
    ad_l = []
    dp_l = []
    oth_l = []
    cells_kept = dict()
    count = 0
    total_cells = 0
    for curr_ind, i in enumerate(indirs):
        print('dir', i)
        ad_f = join(i, f"{prefix_tags_in}.AD.mtx")
        dp_f = join(i, f"{prefix_tags_in}.DP.mtx")
        oth_f = join(i,f"{prefix_tags_in}.OTH.mtx")

        ad, ad_h = load_mtx_df(ad_f, give_header=True)
        dp, dp_h = load_mtx_df(dp_f, give_header=True)
        oth, oth_h = load_mtx_df(oth_f, give_header=True)

        # curr_cell_ids = np.sort(np.array(list(set(dp["Cell"].values).union(set(ad["Cell"].values)).union(set(oth["Cell"].values)))))
        curr_cell_ids = np.arange(1,
                                  max(max(ad_h["Cell"], dp_h["Cell"]),
                                      oth_h["Cell"]) + 1)
        print('curr_cell_ids', curr_cell_ids)
        curr_num_cells = len(curr_cell_ids)
        if num_cells_total is None:
            if cell_subsample != None and (cell_subsample > 0):
                # Subsample a fraction of the cells
                print('curr_num_cells', curr_num_cells)
                num_to_keep = int(
                    round(cell_subsample * curr_num_cells))
                subs = random.choice(curr_cell_ids, replace=False,
                                     size=num_to_keep)
            else:
                raise ValueError(
                    "Either the num_cells_total or cell_subsample need to be filled")
        else:
            # Split the number of cells based on num_cells_total and fraction
            print('is_proportional', is_proportional)
            if is_proportional:
                print('here prop', cell_proportion[curr_ind])
                cells_per_sample = int(
                    round(cell_proportion[curr_ind] * num_cells_total))
            else:
                print('non prop', num_cells_total, len(indirs))
                # Divide evenly
                cells_per_sample = int(
                    round(num_cells_total / len(indirs)))
            if cells_per_sample > curr_num_cells:
                logging.warning("Number of cells less than the desired number. "
                      "Please consider changing a parameter. For now, "
                      "using all cells in the sample.")
                subs = curr_cell_ids
            else:
                subs = random.choice(curr_cell_ids, replace=False,
                                     size=cells_per_sample)

        # Filter
        ad_filt = ad[ad["Cell"].isin(subs)].copy()
        dp_filt = dp[dp["Cell"].isin(subs)].copy()
        oth_filt = oth[oth["Cell"].isin(subs)].copy()
        count += 1

        ##
        # Create the cell coordinates old-new map.
        ##
        # Change the cell coordinates to reflect the #cells in
        # across all samples. e.g. if 2 samples, 100  of each are used,
        # then the first cell in sample 2 will have int coordinate 101
        subs = np.sort(subs)
        cell_coords_map = {val: (ind + total_cells + 1) for ind, val in
                           enumerate(subs)}
        total_cells += len(subs)
        print('subs', subs)
        print('cell coords')
        print(cell_coords_map)
        print('ad_filt')
        print(ad_filt["Cell"])
        ad_filt["Cell"] = ad_filt["Cell"].map(cell_coords_map)
        dp_filt["Cell"] = dp_filt["Cell"].map(cell_coords_map)
        oth_filt["Cell"] = oth_filt["Cell"].map(cell_coords_map)
        print('ad_filt after map')
        print(ad_filt["Cell"])

        # Change the variant coordinates as well after
        ad_filt["Variant"] = ad_filt["Variant"].map(
            vars_coords[i].to_dict())
        dp_filt["Variant"] = dp_filt["Variant"].map(
            vars_coords[i].to_dict())
        oth_filt["Variant"] = oth_filt["Variant"].map(
            vars_coords[i].to_dict())


        cells_kept[i] = cell_coords_map
        ad_l.append(ad_filt)
        dp_l.append(dp_filt)
        oth_l.append(oth_filt)


    # # Since these are filtered, we have to remap again in case
    # # Some variants were removed.
    # new_vars = set()
    # for ind, val in enumerate(indirs):
    #     new_vars = new_vars.union(set(ad_l[ind]["Variant"].values))
    #     #new_vars.union(set(vars_coords[ind]))
    #
    # # Sort the new variants
    # new_vars = np.sort(np.array(list(new_vars)))
    # # Map the new vars and update the values
    # new_vars_map = {x:ind+1 for ind, x in enumerate(new_vars)}
    # for ind, val in enumerate(indirs):
    #     ad_l[ind]["Variant"] = ad_l[ind]["Variant"].map(new_vars_map).astype(int)
    #     dp_l[ind]["Variant"] = dp_l[ind]["Variant"].map(new_vars_map).astype(int)
    #     oth_l[ind]["Variant"] = oth_l[ind]["Variant"].map(new_vars_map).astype(int)

    ####
    # 5. Merge the dataframe matrices
    ####
    ad_full = pd.concat(ad_l, axis=0)
    print('ad_full')
    print(ad_full.head())
    dp_full = pd.concat(dp_l)
    oth_full = pd.concat(oth_l)

    ####
    # . Save the matrices and the labels
    ####
    # Save the order of IDs and the cell maps
    # Save the pseudo matrices
    # Need to also add in a row tha is max var, max cell, number of entries
    wrap_write_mtx_df(outdir, ad_full, dp_full, oth=oth_full, to_rm=True,
                      prefix=prefix_tags_out)

    # Save cell indices
    for ind, val in enumerate(indirs):
        with open(join(outdir, f"cell_indices_{ind}.txt"), 'w') as f:
            curr = val + "\n" + "old index,new index"
            for k in cells_kept[val]:
                curr = f"{curr}\n{k},{cells_kept[val][k]}"
            f.write(curr)

    in_files = [os.path.join(val, f"{prefix_samples_in}.samples.tsv") for val in indirs]
    logging.info('in_files')
    logging.info(in_files)
    add_suffix_to_labels(in_files, cells_kept,
                         join(outdir, "cells_meta.tsv"),
                         sample_names=sample_names)
    return


def merge_vcf(vcf_files, out_vcf):
    cmd = f"""bcftools merge {vcf_files} - Oz - o {out_vcf}"""
    print(cmd)
    os.system(cmd)
    return



# Subsampling on bam files:
# Pros: More raw
# Cons: Will subsampling on bam miss the paired reads? Can that be connected

#Subsampling on mtx files:
#Pros: Already takes care of the bam processing, and takes care of batch
# for each separately.
@click.command()
@click.argument("outdir", type=click.Path())
@click.argument("indirs", type=click.Path(exists=True), nargs=-1)
@click.option("--is_prop", default=False, type=click.BOOL)
@click.option("--num_cells", default=1000)
@click.option("--samples", default="")
@click.option("--prefix_in", default="")
@click.option("--prefix_out", default="")
@click.option("--prefix_tags_in", default="cellSNP.tag")
@click.option("--prefix_vcf_in", default="cellSNP.base")
@click.option("--prefix_tags_out", default="cellSNP.tag")
@click.option("--prefix_vcf_out", default="cellSNP.base")
@click.option("--prefix_samples_in", default="cellSNP")
def main(outdir, indirs, is_prop, num_cells, samples,prefix_in,prefix_out,
         prefix_tags_in, prefix_vcf_in, prefix_tags_out, prefix_vcf_out,
         prefix_samples_in):
    logging.basicConfig(level=logging.DEBUG)
    logging.info('samples')
    logging.info(samples)
    logging.info(type(samples))
    if samples != "":
        samples = samples.split(',')
    if samples == []:
        samples=None
    logging.info(type(samples))

    # Precedence of prefix if prefix_in and/or prefix_out are set
    if prefix_in != "":
        prefix_tags_in = prefix_in,
        prefix_vcf_in = prefix_in
        prefix_samples_in = prefix_in
    if prefix_out != "":
        prefix_tags_out = prefix_out
        prefix_vcf_out = prefix_out


    subsample_sparse_matrices(outdir, indirs,
                              num_cells_total=num_cells,
                              is_proportional=is_prop,
                              cell_subsample=-1, sample_names=samples,
                              prefix_tags_in=prefix_tags_in,
                              prefix_vcf_in=prefix_vcf_in,
                              prefix_tags_out=prefix_tags_out,
                              prefix_vcf_out=prefix_vcf_out,
                              prefix_samples_in=prefix_samples_in)
    return


if __name__ == "__main__":
    main()
#     os.chdir("/data2/mito_lineage/Analysis/multiplex")
#     main(["data/vireo/pseudo/numC1000_ispropFalse",
#          "data/PBMC_J_cellSNP", "data/PBMC_P_cellSNP",
#          '--num_cells', 1000, '--is_prop', False])
