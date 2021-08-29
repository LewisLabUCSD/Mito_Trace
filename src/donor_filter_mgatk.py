""" Snakemake Script. Don't import snakemake."""

import pandas as pd
from os.path import dirname, join
import click

print('input', snakemake.input)
print('output ', snakemake.output)
#for sample, infile in zip(snakemake.params.samples, snakemake.input):
cells_meta = pd.read_csv(snakemake.params.cells_meta , sep='\t')
sample = snakemake.params.sample

for d, df in cells_meta.groupby('donor'):
    d = int(d)
    nt_pileups = {}
    cond_positions = {}
    nt_pileups["coverage"] = pd.DataFrame()

    for ind, curr_cov_f in enumerate(snakemake.input.cov):
        print('curr_cov_f ', curr_cov_f)
        curr_samp_cov = pd.read_csv(curr_cov_f,
                                    header=None)  # load sample coverage
        condition = sample[ind]
        curr_samp_donor_df = df[
            df["condition"] == condition]  # Condition-donors
        filt_df = curr_samp_cov.loc[
            curr_samp_cov[1].isin(curr_samp_donor_df["raw ID"])]
        # print('filt_df', filt_df.head())
        filt_df[1] = filt_df[1] + "_" + sample[ind]
        cond_positions[sample[ind]] = set(filt_df[0].values)
        nt_pileups["coverage"] = nt_pileups["coverage"].append(filt_df,
                                                               ignore_index=True)

    # Filter by overlapping positions:
    pos_to_keep = list(cond_positions.values())[0].intersection(
        *list(cond_positions.values()))
    nt_pileups["coverage"] = nt_pileups["coverage"][
        nt_pileups["coverage"][0].isin(pos_to_keep)]

    # Some cells may not be present anymore, need to store them.
    cells_to_keep = set(nt_pileups["coverage"][1])

    # Save coverage and depth
    nt_pileups["coverage"].to_csv(str(snakemake.output[int(d)]), index=False,
                                  header=False)
    depth = nt_pileups["coverage"].groupby(1).sum()[2] / 16569
    depth.to_csv(join(dirname(snakemake.output[int(d)]), f"d{d}.depthTable.txt"),
                 sep='\t', header=False)

    # Filter for the nucleotides as well
    for nuc in ["C", "A", "G", "T"]:
        nt_pileups[nuc] = pd.DataFrame()
        for ind, curr_cov_f in enumerate(snakemake.input.cov):
            curr_samp_cov = pd.read_csv(join(dirname(curr_cov_f),
                                             f"{sample[ind]}.{nuc}.txt"),
                                        header=None)  # load sample coverage
            condition = sample[ind]
            curr_samp_donor_df = df[df["condition"] == condition]
            filt_df = curr_samp_cov.loc[curr_samp_cov[1].isin(
                curr_samp_donor_df["raw ID"])]  # Filter by donor
            ## Filter by positions to keep and cells
            filt_df[1] = filt_df[1] + "_" + sample[ind]
            filt_df = filt_df[((filt_df[0].isin(pos_to_keep)) & (
                filt_df[1].isin(cells_to_keep)))]

            nt_pileups[nuc] = nt_pileups[nuc].append(filt_df,
                                                     ignore_index=True)
        # Save coverage
        nt_pileups[nuc].to_csv(
            join(dirname(snakemake.output[int(d)]), f"d{d}.{nuc}.txt"),
            header=False, index=False)

    cells_meta = cells_meta[cells_meta["ID"].isin(cells_to_keep)]
    cells_meta.to_csv(join(dirname(snakemake.output[int(d)]), "cells_meta.tsv"),
                      sep='\t', index=False)

