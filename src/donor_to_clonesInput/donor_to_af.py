""" Snakemake Script. Don't need to import snakemake."""

import pandas as pd
from os.path import dirname, join
import numpy as np
#import click
#import snakemake
from src.utils.data_io import af_to_vireo

print('input', snakemake.input)
print('output ', snakemake.output)
print('params', snakemake.params)

cells_meta = pd.read_csv(snakemake.params.cells_meta , sep='\t')
sample = snakemake.params.sample
ref_fa = pd.read_csv(snakemake.params.ref_mt, sep='\t',
                     header=None, index_col=0)


for d, curr_donor in cells_meta.groupby('donor'):
    d = int(d)
    print('donor', d)
    print(dirname(snakemake.output[d]))
    vars_conds = []
    af_donor_conds = []
    dp_donor_conds = []

    # Loop through each condition and filter cells
    for ind, curr_cov_f in enumerate(snakemake.input.cov):
        print('curr_cov_f ', curr_cov_f)
        ## AF
        af = pd.read_csv(join(dirname(curr_cov_f), "af_by_cell.tsv"), sep='\t',
                         header=0, index_col=0) #cells-by-variants
        condition = sample[ind]
        af.index = af.index + "_" + condition
        #print('af')
        #print(af.head())
        curr_samp_donor_df = curr_donor[
            curr_donor["condition"] == condition]  # Condition-donors
        af_donor_samp = af.loc[af.index.isin(curr_samp_donor_df["ID"].values)]
        vars_conds.append(set(af_donor_samp.columns))
        af_donor_conds.append(af_donor_samp)

        ## Coverage
        curr_samp_cov = pd.read_csv(curr_cov_f,
                                    header=None)  # load sample coverage
        curr_samp_cov[1] = curr_samp_cov[1] + "_" + sample[ind]
        curr_samp_cov = curr_samp_cov.loc[
            curr_samp_cov[1].isin(curr_samp_donor_df["ID"])]
        curr_samp_cov = curr_samp_cov.rename({1: "Cell", 0: "Position", 2: "Value"}, axis=1).pivot(index="Cell",
                                                                                                  columns="Position",
                                                                                                  values="Value")
        dp_donor_samp = pd.DataFrame(index=af_donor_samp.index, columns=af_donor_samp.columns)

        def fill_dp(x,y):
            #print('x', x)
            #print('y', y.head())
            return y.loc[y.index.isin(x.index), int(x.name[:-1])].values

        dp_donor_samp = dp_donor_samp.apply(fill_dp , args=(curr_samp_cov,))

        dp_donor_samp.isnull()
        #dp_donor_samp = dp_donor_samp.loc[dp_donor_samp.index.isin(af_donor_samp.index), dp_donor_samp.columns.isin([int(x[:-1] ) for x in af_donor_samp.columns])]
        dp_donor_conds.append(dp_donor_samp)

    # Get overlap variants and remove cells/vars with 0s
    overlap_vars = vars_conds[0].intersection(*vars_conds)
    af_donor_df = pd.concat(af_donor_conds, join='inner', axis=0,
                                                ignore_index=False)
    af_donor_df = af_donor_df.loc[(af_donor_df.sum(axis=1)>0), (af_donor_df.sum(axis=0)>0)]
    #print('af_donor_df')
    #print(af_donor_df.head())
    dp_donor_df = pd.concat(dp_donor_conds, join='inner', axis=0,
                                                ignore_index=False).fillna(0)
    dp_donor_df = dp_donor_df.loc[dp_donor_df.index.isin(af_donor_df.index), dp_donor_df.columns.isin(af_donor_df.columns)]
    #print('dp_donor_df')
    #print(dp_donor_df.head())
    print("af null")
    print(af_donor_df.isnull().sum().sum())
    print("dp null")
    print(dp_donor_df.isnull().sum().sum())

    print(dp_donor_df.shape)
    print(af_donor_df.shape)
    # Save allele + depth matrix + cells_meta
    af_donor_df.to_csv(join(dirname(snakemake.output[int(d)]), "af.tsv"),
                      sep='\t', index=True)
    dp_donor_df.to_csv(join(dirname(snakemake.output[int(d)]), "dp.tsv"),
                      sep='\t', index=True)

    curr_cells_meta = cells_meta[cells_meta["ID"].isin(af_donor_df.index)].copy()
    curr_cells_meta.to_csv(join(dirname(snakemake.output[int(d)]), "cells_meta.tsv"),
                      sep='\t', index=False)

    # Create VireoIn
    # Requires AD sparse matrix
    af_to_vireo(af_donor_df.transpose(), dp_donor_df.transpose(),
                outdir=dirname(snakemake.output[int(d)]),
                out_name="cellSNP")

    curr_vars = af_donor_df.columns
    import pandas as pd

    af_meta = pd.DataFrame({ "POS": [int(x[:-1]) for x in curr_vars],
                                                    "ALT": [x[-1] for x in curr_vars],
                                                    "index": np.arange(1,len(curr_vars)+1),
                                                    "REF": [f"{ref_fa.loc[int(x[:-1]),1]}_{x[:-1]}" for x in curr_vars]
                                                }
                                            )
    af_meta["#CHROM"] = "MT"
    af_meta[["#CHROM", "POS", "REF","ALT", "index"]].to_csv(join(dirname(snakemake.output[d]), "cellSNP.base.vcf"),
                   sep="\t", index=False)
