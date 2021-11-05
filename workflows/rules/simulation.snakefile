wildcard_constraints:
    num="\d+"

import os
import pandas as pd
from snakemake.utils import validate



rule all:
    input:
        expand("figures/performance/{num}_lineage.png",num=config['num']),
        expand("figures/performance/{num}_growth.png",num=config['num']),
        "figures/performance/all_growth_compare.png"
        #expand("figures/performance/{num}_enrichment.png",num=config['num'])
# -----------------------------------

rule simulate_initial:
    #"{num}_params.yaml"
    output:
        init_cell_af = "data/simulate/{num}_cell_af.csv",
        init_cell_clone = "data/simulate/{num}_cell_clone.csv",
        MT_clone_map = "data/simulate/{num}_mt_clone.csv"

rule simulate_sequence:
    input:
        init_cell_af = "data/simulate/{num}_cell_af.csv",
        init_cell_clone = "data/simulate/{num}_cell_clone.csv",
    output:
        init_cell_seq = "data/simulate/{num}_cell_af_seq.csv",


rule simulate_growth:
    input:
        init_cell_af = "data/simulate/{num}_cell_af.csv",
        init_cell_clone = "data/simulate/{num}_cell_clone.csv",
        MT_clone_map = "data/simulate/{num}_mt_clone.csv"
    params: type = config['growth']
    output:
        grow_cell_af = temp("data/simulate/{num}_grow_cell_af.csv"),
        grow_cell_clone = temp("data/simulate/{num}_grow_cell_clone.csv"),
        grow_MT_clone_map = "data/simulate/{num}_mt_clone.csv"
    #shell:

rule simulate_growth_sequence:
    input:
        grow_cell_af = "data/simulate/{num}_grow_cell_af.csv",
        grow_cell_clone = "data/simulate/{num}_grow_cell_clone.csv"
    params: n_subsample = config["n_subsample"]
    output:
        subsample_cell_af = "data/simulate/{num}_subs_cell_af.csv",
        subs_cell_clone_meta = "data/simulate/{num}_subs_cell_clone.csv"


rule model_call_variants_threshold:
    """ Call variants using one of a few types. 
    
    :type type = {'thresh', 'GMM', 'mgatk'}
    Modes: 
        thresh: Sets a simple threshold for number of cells and reads a
                variant needs to be called a variant. 
        GMM: Gaussian mixture model for quality of the variants. 
        mgatk: This will use strand concordance, allele quality, and 
            number of reads to determine the correct threshold. This only
            works information on aligned read strands in scATAC. 
         
    Include a model for known variants, such as transversion"""
    input:
        subsample_cell_af = "data/simulate/{num}_subs_cell_af.csv",
    params: type = config['variants']['type']
    output:
        filt_cell_af = "data/model/{num}_varfilt_cell_af.csv", # If making it binary, could be better as a list instead of matrix for memory
        variants = "data/model/{num}_pred_vars.csv",


rule model_call_lineages:
    input:
        rules.model_call_variants_threshold.output.filt_cell_af
    output:
        lineage_pred = "data/model/{num}_lineage_labels.csv"


rule performance_lineages:
    """ Looks at performance for cell-clone relationship
    and MTvariant-clone relationship
    
    If a tree and not cluster labels are given, will assume MT lineages
    are based on variants, 
    """
    input:
        labels = rules.model_call_lineages.output.lineage_pred,
        true = rules.simulate_growth_sequence.output.subs_cell_clone_meta
    output:
        lin_perf = "figures/performance/{num}_lineage.png",
        var_perf = "results/lineage/{num}_variantcall_performance.png"


rule model_lineage_growth_estimate:
    """ Will estimate the growth rate based on the number of cells in
    the before and after experiments.
    
    """
    input:
        rules.simulate_growth_sequence.output.subs_cell_clone_meta,
        rules.model_call_lineages.output
    output:
        grow_vals = "results/growth/{num}_growth_estimate.csv",
        grow_fig = "results/growth/{num}_growth_estimate.png"



rule performance_growth_estimate:
    input:
        rules.model_lineage_growth_estimate.output.grow_vals
    params: config['growth']
    output:
        "figures/performance/{num}_growth.png"


rule model_agg_growth_phenotype_model:
    """ Compare the growth of the largest clones between the disease and
    non-disease samples. 
    
    """
    input: expand("results/growth/{num}_growth_estimate.csv", num=config['num'])
    params: config["samples"] # Which sample is which value
    output:
        samples_fig = "figures/performance/all_growth_compare.png",
        samples_pred = "results/compare_growth/growth_distribution.csv"


# rule model_lineage_enrichment:
#     input:
#
#     output:
#     shell:



# rule performance_enrichment:
#     input:
#     output:

