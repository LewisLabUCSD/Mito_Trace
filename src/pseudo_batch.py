from os.path import join
import os
from numpy import random
import numpy as np
import pandas as pd
import click
import sys


def run_batch(indirs, outdir, num_reads_total):
    num_reads_each = len(indirs)/num_reads_total
    return


def load_mtx_df(in_f, skip_first=True, give_header=False):
    df = pd.read_csv(in_f, comment="%", header=None, sep="\t")
    df.columns = ["Variant", "Cell", "integer"]
    if skip_first:
        head = df.iloc[0]
        df = df.iloc[1:] # Seems to be summary values
    if give_header:
        return df,head
    return df


def merge_vcf_ids(indirs, outdir=None):
    """ Merges the variant files called from cellSNP

    :param indirs:
    :return:
    """

    curr_vcf = join(indirs[0], "cellSNP.base.vcf")
    if not os.path.exists(curr_vcf):
        if not os.path.exists(curr_vcf + ".gz"):
            raise ValueError("VCF file not here!")
        else:
            curr_vcf = curr_vcf + ".gz"
    variants = pd.read_csv(curr_vcf, skiprows=1, sep="\t")
    variants["old " + indirs[0]] = variants.index.values.astype(int) + 1
    print('variants')
    print(variants.head())
    print(variants.tail())
    old_variants = {}
    old_variants[indirs[0]] = variants.copy()
    old_variants[indirs[0]].loc[:, "old"] = 0
    old_variants[indirs[0]].loc[:, "old"] = old_variants[
                                                indirs[0]].index + 1
    # variants[indirs[0]] = variants.index.values + 1
    for ind, val in enumerate(indirs[1:]):
        curr_vcf = join(val, "cellSNP.base.vcf")
        if not os.path.exists(curr_vcf):
            if not os.path.exists(curr_vcf + ".gz"):
                raise ValueError("VCF file not here!")
            else:
                curr_vcf = curr_vcf + ".gz"
        curr = pd.read_csv(curr_vcf, skiprows=1,
                           sep="\t")
        # curr["old index"] = curr.index.values+1
        # curr[val] = curr.index.values + 1
        curr["old " + val] = curr.index.values.astype(int) + 1
        variants = pd.merge(variants, curr, on=["#CHROM", "POS", "ALT"],
                            how="outer")

        print('variants')
        print(variants.head())
        print(variants.tail())
        old_variants[val] = curr.copy()
        old_variants[val]["old"] = 0
        old_variants[val]["old"] = curr.index.values + 1

    # Loop again and map the coordinates
    vars_coords = dict()
    variants["new ID"] = variants.index.values + 1
    full_vars = variants.copy()

    # what the new index should be.
    for val in indirs:
        curr_vcf = join(val, "cellSNP.base.vcf")
        if not os.path.exists(curr_vcf):
            if not os.path.exists(curr_vcf + ".gz"):
                raise ValueError("VCF file not here!")
            else:
                curr_vcf = curr_vcf + ".gz"
        curr_vars = pd.read_csv(curr_vcf, skiprows=1, sep="\t")
        curr_vars["old index"] = curr_vars.index + 1
        curr_vars = pd.merge(curr_vars, full_vars, how="inner",
                             on=["#CHROM", "POS", "ALT"])
        curr_vars.index = curr_vars["old index"]
        vars_coords[val] = curr_vars["new ID"]  # Pandas series

    if outdir is not None:
        variants.to_csv(join(outdir, "cellSNP.base.vcf"), sep='\t',
                        index=False)
        for ind, val in enumerate(indirs):
            out_f = join(outdir, f"variant_indices_{ind}.tsv")
            vars_coords[val].index.name = val
            vars_coords[val].to_csv(out_f,
                                    sep="\t")  # Use index + header
    return vars_coords, variants



def subsample_sparse_matrices(outdir, indirs, cell_subsample=0.1,
                              num_cells_total=1000,
                              is_proportional=False):
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
    # Count number of cells first
    cell_count = {}
    for i in indirs:
        curr_dp_f = join(i, "cellSNP.tag.DP.mtx")
        curr_ad_f = join(i, "cellSNP.tag.AD.mtx")
        curr_oth_f = join(i, "cellSNP.tag.OTH.mtx")
        dp = load_mtx_df(curr_dp_f)
        ad = load_mtx_df(curr_ad_f)
        oth = load_mtx_df(curr_oth_f)
        if len(ad) == 0:
            raise ValueError("The allele matrix is empty. It could have all been filtered out.")
        # cell_count[i] = len(set(dp["Cell"].values))
        cell_count[i] = len(list(
            set(dp["Cell"].values).union(set(ad["Cell"].values)).union(
                set(oth["Cell"].values))))
    cell_proportion = np.array(list(cell_count.values())) / sum(
        cell_count.values())
    print('cell proportion', cell_proportion)

    vars_coords, variants = merge_vcf_ids(indirs, outdir=outdir)
    # ad["Variant"] = vars_coords[ad["Variant"]]
    # dp["Variant"] = vars_coords[dp["Variant"]]
    # oth["Variant"] = vars_coords[oth["Variant"]]

    ad_l = []
    dp_l = []
    oth_l = []
    cells_kept = dict()
    count = 0
    total_cells = 0
    for i in indirs:
        print('dir', i)
        ad_f = join(i, "cellSNP.tag.AD.mtx")
        dp_f = join(i, "cellSNP.tag.DP.mtx")
        oth_f = join(i, "cellSNP.tag.OTH.mtx")

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
                print('here prop', cell_proportion[count])
                cells_per_sample = int(
                    round(cell_proportion[count] * num_cells_total))
            else:
                print('non prop', num_cells_total, len(indirs))
                # Divide evenly
                cells_per_sample = int(
                    round(num_cells_total / len(indirs)))
            if cells_per_sample > curr_num_cells:
                print("Number of cells less than the desired number. "
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
        # Merge the samples together
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

        # merging with merge_vcf_ids above
        # subsample_d["cell"].append(set(cell_coords_map.values()))

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

    # For each list of dataframes, merge and save
    ad_full = pd.concat(ad_l, axis=0)
    print('ad_full')
    print(ad_full.head())
    dp_full = pd.concat(dp_l)
    oth_full = pd.concat(oth_l)

    # Save the order of IDs and the cell maps

    # Save the pseudo matrices
    # Need to also add in a row tha is max var, max cell, number of entries
    header = "%%MatrixMarket matrix coordinate integer general\n%\n"

    if os.path.exists(join(outdir, "cellSNP.tag.DP.mtx")):
        os.remove(join(outdir, "cellSNP.tag.DP.mtx"))
    if os.path.exists(join(outdir, "cellSNP.tag.AD.mtx")):
        os.remove(join(outdir, "cellSNP.tag.AD.mtx"))
    if os.path.exists(join(outdir, "cellSNP.tag.OTH.mtx")):
        os.remove(join(outdir, "cellSNP.tag.OTH.mtx"))

    with open(join(outdir, "cellSNP.tag.AD.mtx"), 'a') as file:
        file.write(header)
        ad_full = pd.concat((pd.DataFrame(
            {"Variant": ad_full["Variant"].max(),
             "Cell": ad_full["Cell"].max(),
             "integer": ad_full.shape[0]}, index=["Meta"]),
                             ad_full.sort_values(["Variant", "Cell"])))
        ad_full.to_csv(file, sep="\t", header=False, index=False)

    with open(join(outdir, "cellSNP.tag.DP.mtx"), 'a') as file:
        file.write(header)
        dp_full = pd.concat((pd.DataFrame(
            {"Variant": dp_full["Variant"].max(),
             "Cell": dp_full["Cell"].max(),
             "integer": dp_full.shape[0]}, index=["Meta"]),
                             dp_full.sort_values(["Variant", "Cell"])))
        dp_full.to_csv(file, sep="\t", header=False, index=False)
    with open(join(outdir, "cellSNP.tag.OTH.mtx"), 'a') as file:
        file.write(header)
        oth_full = pd.concat((pd.DataFrame(
            {"Variant": oth_full["Variant"].max(),
             "Cell": oth_full["Cell"].max(),
             "integer": oth_full.shape[0]}, index=["Meta"]),
                              oth_full.sort_values(
                                  ["Variant", "Cell"])))
        oth_full.to_csv(file, sep="\t", header=False, index=False)

    # Save cell indices
    for ind, val in enumerate(indirs):
        with open(join(outdir, f"cell_indices_{ind}.txt"), 'w') as f:
            curr = val + "\n" + "old index,new index"
            for k in cells_kept[val]:
                curr = f"{curr}\n{k},{cells_kept[val][k]}"
            f.write(curr)
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
def main(outdir, indirs, is_prop, num_cells):
    subsample_sparse_matrices(outdir, indirs,
                              num_cells_total=num_cells,
                              is_proportional=is_prop,
                              cell_subsample=0.1)
    return


if __name__ == "__main__":
    main()
#     os.chdir("/data2/mito_lineage/Analysis/multiplex")
#     main(["data/vireo/pseudo/numC1000_ispropFalse",
#          "data/PBMC_J_cellSNP", "data/PBMC_P_cellSNP",
#          '--num_cells', 1000, '--is_prop', False])
