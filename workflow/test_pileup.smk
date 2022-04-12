report: "report/workflow.rst"

#from src.config import ROOT_DIR
#workdir: ROOT_DIR
wildcard_constraints:
    cellr='True|False',
    kparam='[0-9]+',  #'%d', #[0-9]+',
    variants = "simple|mgatkdonor|init|simpleUnion",
    d = "[0-9]+"

#configfile:  "parameters/pipeline/cosmo_server/jan21_2021.yaml"
from src.config import ROOT_DIR
from src.utils.parse_config import read_config_file
import os
import numpy as np
from os.path import join, dirname
import pandas as pd
from snakemake.utils import min_version
from icecream import ic
min_version("6.0")
print('config', config)

import pickle
from snakemake.utils import Paramspace

cells_pspace = Paramspace(config["cells_setup"]["params"], filename_params="*", param_sep="__")
seq_pspace = Paramspace(config["seq_setup"]["params"], param_sep="__")

########################################################################
# Setup parameters and outdir
########################################################################
res = join(config["outdir"], "clone_pileups_simulation", config["prefix"])

params = read_config_file(config["config"])
samples = pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
#anno_res = join(config["outdir"], "annotation", "data", params['annotations']['version'], config["prefix"])
anno_res = join(config["outdir"], "annotation", "data", config["prefix"])


# Merge the experiment-specific config with the pipeline parameters
params["results"] = res
#params["samples_meta"] = config["samples_meta"]
params["samples"] = samples
params["genome"] = config["genome"]
params["N_DONORS"] = config["N_DONORS"]
params["anno_res"] = anno_res
config["samples"] = samples
config["params"] = params
####

params_mt = params["mtpreproc"]
params_ft = params["filters"]
ft = params_ft["params"]
params_mult = params["multiplex"]
params_mgatk = params["mgatk"]
params_clones = params["clones"]
params_annclo = params["annotation_clones"]["params"]
params_clch = params["de_clones_change"]["params"]
params_clch["samples"] = samples

cellrbc = params_mt["params"]["cellrbc"]
num_reads_filter = params_mt["params"]["numreadsfilter"]
maxBP = params_mt["maxBP"]
ref_fa = params["genome_path"][config["genome"]]["ref_fa"]
mt_ref_fa = params["genome_path"][config["genome"]]["mt_ref_fa"]

gff = params["genome_path"][config["genome"]]["gff"]


ncellsthreshmgatk = params_mgatk["params"]["ncellsthresh"]
num_cells = params_mult["pseudo_multiplex"]["num_cells"]
is_prop = params_mult["pseudo_multiplex"]["is_proportional"]

nclonelist = params_clones['vireo']['params']['nclonelist']
config["params_clones"] = params_clones



###################################
## Rules and workflow
###################################
###################################
rule all:
    input:
        expand("{outdir}/{sim_cells}/{sim_seq}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf",
            out_dir=res,sim_cells=cells_pspace.instance_patterns, sim_seq=seq_pspace.instance_patterns,
            cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),

        # Multiplex: Multiplexed VAF and depth dendrograms for each donor and selected variants
        expand("{out_dir}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),

        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            method=params_clones["method"],
            kparam=params_clones["knn"]["params"]["resolution"],d=np.arange(config["N_DONORS"])),

        # Variants
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/anno_variants.tsv",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh']),

        # Clone stats
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/_clone_complete.txt",
            out_dir=res, cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            method=params_clones["method"],
            kparam=params_clones["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x!="simple"]
        ),

        # Enrichment
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/_enrichment_complete",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            variants=[x for x in params_clones["variants"] if x != "simple"],
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
            kparam=params_clones["knn"]["params"]["resolution"]),

        # ##  Methods compare
        # expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/comparisons.ipynb",
        #     out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
        #     mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
        #     hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),
        #

################################################################
## Import from prior snakefile modules
## Here, we redefine the input to be based on our config['files'] dictionary
from snakemake.utils import min_version
min_version("6.0")
module mtpreprocMod:
    snakefile: "./rules/mt_preprocess.smk"
    config: params

module mgatkMod:
    snakefile: "./rules/mgatk.smk"
    config: params

module multMod:
    snakefile: "./rules/multiplex.smk"
    config: params



########################################################################
## 1. Go from 10x output to filtered scPileup outdir.
########################################################################
# rule filter_cell_bc:
#     """Extracts only the relevant cell barcodes and removes the .bam from the barcode names."""
#     input:
#         all = expand("{{outdir}}/data/{{sample}}/MT/scPileup_concat_{{num_read}}/numread_{{num_read}}_all.{nt}.strands.txt",
#                      nt=["coverage", "A", "C", "G", "T"]),
#         barcode_p = "{outdir}/data/{sample}/MT/cellr_{cellrbc}/{sample}_barcode_data.p"
#     output:
#         "{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt"
#     run:
#         for n in ["A", "C", "G", "T", "coverage"]:
#             #print('n', n)
#             curr_f = input.all[0]
#             curr_f = curr_f.replace(".coverage.", "." + n + ".")
#             #print('curr_f', curr_f)
#             df = pd.read_csv(curr_f, sep=',')
#             df["CB"] = df["CB"].str.replace(".bam","")
#             barcodes = pickle.load(open(input.barcode_p, "rb"))
#             if isinstance(barcodes, dict):
#                 df = df[df["CB"].isin(list(barcodes.keys()))]
#             else:
#                 df = df[df["CB"].isin(barcodes)]
#             curr_out_f = output[0]
#             curr_out_f = curr_out_f.replace(".coverage.", "." + n + ".")
#             if n == "coverage":
#                 df = df.iloc[:,:3]
#             df = df.sort_values(["CB", "Position"])
#             df.to_csv(curr_out_f, header=None, index=None)
#
#
# rule scPileup_MT_matrix:
#     """Create the position-by-cell coverage matrix"""
#     input:
#         all = rules.scPileup_mergeStrands.output.all[0],  #"{outdir}/data/{sample}/MT/scPileup_concat_{num_read}/numread_{num_read}_all.coverage.strands.txt.gz",
#         barcode_p = rules.barcode_filter.output[0] #"{outdir}/data/{sample}/MT/cellr_{cellrbc}/{sample}_barcode_data.p"
#     output:
#         sc_coverage_f = "{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/sc_coverage.csv"
#     threads: 16
#     params:
#         maxBP= 16571
#     shell:
#         "python src/figures/plot_heatmap_coverage.py sc_mt {input.barcode_p} {input.all} {output.sc_coverage_f} {params.maxBP}"
#
#
# rule plot_sc_coverageBar_and_heat:
#     """Plot the posiitonal coverages."""
#     input:
#         sc_coverage_f = "{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/sc_coverage.csv"
#     output:
#         save_f_coverage = report("{outdir}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png", category="mtpreproc"),
#         save_f_heat = report("{outdir}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position.png", category="mtpreproc"),
#     shell:
#         "python src/figures/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_heat} {output.save_f_coverage}"


def get_filt(w):
    return w.mincells, w.minreads, w.topN, w.hetthresh, w.minhetcells, w.hetcountthresh, w.bqthresh


rule create_filters:
    input:
        concat_dir =  "{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt" #(rules.filter_cell_bc.output[0]) #"{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt.gz"
    output:
        cov = "{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/{sample}.coverage.txt",
        af = "{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/af_by_cell.tsv",
        fig = report("{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/heatmap.png",
                      category="mtpreproc", subcategory="filter"),
        fig2 = report("{outdir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/initial_cell_depth.png",
                      category="mtpreproc",subcategory="filter"),
    params:
        concat_d = lambda wildcards, input: dirname(input.concat_dir),
        ref_fa = params["genome_path"][config['genome']]['mt_ref_fa'],
        name = lambda wildcards: wildcards.sample,
        filt_params = get_filt,
    resources:
        mem_mb=90000
    shell: "python src/mtpreproc/calculate_AF_by_cell.py {params.concat_d} {output.af} {params.ref_fa} {params.name} {params.filt_params}" # --log {log}"


########################################################################
## 2. Call variants using MGATK and convert to the Vireo multiplex input
########################################################################
use rule * from mgatkMod

use rule mgatk_to_vireoIn from mgatkMod with:
    input:
        "{out_dir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/{sample}.variant.rds"
    output:
        "{out_dir}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx"


rule mgatk_to_variant:
    input:
        cells = lambda wildcards: expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx",
                                  sample=samples.index)
    output:
        vcf = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/variants.vcf",
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "variant_types", "scoreVariants.ipynb"),
        vcfs = lambda wildcards, input: [join(dirname(x),"cellSNP.base.vcf") for x in input.cells],
        mt_fasta = mt_ref_fa
    run:
        #all_vcf = pd.read_csv(params.vcfs[0], sep="\t")[["#CHROM", "POS", "REF", "ALT"]]
        #for i in params.vcfs[1:]:
            #all_vcf = all_vcf.merge(pd.read_csv(i, sep="\t")[["#CHROM", "POS", "REF", "ALT"]], how="inner")
        all_vcf = []
        for i in params.vcfs:
            all_vcf.append(pd.read_csv(i, sep="\t")[["#CHROM", "POS", "REF", "ALT"]])
        all_vcf = pd.concat(all_vcf, ignore_index=True).drop_duplicates()
        all_vcf["QUAL"] = "."
        all_vcf["FILTER"] = "."
        all_vcf["INFO"] = "."
        all_vcf["ID"] =  all_vcf["REF"] +  ">" + all_vcf["ALT"]
        all_vcf["REF"] = [x[-1] for x in all_vcf["REF"]]
        header = "##fileformat=VCFv4.0"
        header = header + "\n" + f"##reference=file:/{params.mt_fasta}\n"
        print('header', header)
        print('output.vcf', output.vcf)
        with open(output.vcf, 'a') as file:
            file.write(header)
        print(all_vcf.columns)
        all_vcf.to_csv(output.vcf, mode='a', index=False, sep="\t", header=True)

rule merge_variant_samples:
    input:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/variants.vcf",
    output:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/anno_variants.tsv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/anno_variants/anno_variants.ipynb"
    params:
        note = join(ROOT_DIR, "workflow", "notebooks", "variant_types", "scoreVariants.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
    shell:
        "papermill -p vcf {input[0]} -p outdir {params.outdir} {params.note} {output.note}"

########################################################################
## 3. Merge the sample inputs together by concatenating.
########################################################################
##### The cell barcodes get changed here!
rule merged:
    """ Merge the sample pileup matrices
    """
    input:
        lambda wildcards: expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/{sample}/vireoIn/cellSNP.tag.AD.mtx",
                                        sample=samples.index)# , minC=wildcards.minC, minAF=wildcards.minAF)
    output: "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/cellSNP.tag.AD.mtx"
    #input: lambda wildcards: expand("data/{{prefix}}/chrM/{name}_cellSNP_minC{{mt_minC}}_minAF{{mt_minAF}}", name=config["samples"])
    params:
        num_cells=num_cells,
        is_prop = is_prop,
        indirs = lambda wildcards, input: [dirname(x) for x in input],
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix= ','.join(samples.index)
    #log: "logs/{prefix}/vireo/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}.log"
    shell: "python -m src.pseudo_batch {params.outdir} {params.indirs} --num_cells {params.num_cells} --is_prop {params.is_prop} --samples {params.prefix}" # > {log} 2>&1"


########################################################################
## 4. De-multiplex the merged conditions using Vireo
########################################################################
use rule * from multMod as mult_*

use rule donors_type_variants from multMod with:
    output:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/variants/variants.ipynb",


use rule donors_plotAF from multMod as mult_donors_plotAF with:
    input:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/dendrograms/figures/donor{d}_dendrogram.{suf}",
               d=np.arange(config["N_DONORS"]), suf=["png", "depth.png", "withHigh.png"]),
        #"{out_dir}/merged/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/dendrograms/af_dendro.ipynb",

rule multiplex_report:
    input:
        multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/",
                 "multiplex_AF_SNPs_all_afFilt.png", "multiplex_clusters_all.labels.png")
    output:
        report("{out_dir}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/multiplex.pdf", category="Multiplex")
    shell:
        "convert {input} {output[0]}"

rule donor_report:
    input:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/dendrograms/figures/donor{{d}}_dendrogram.{suf}",
               suf=["png", "depth.png", "withHigh.png"]),
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/mgatk/d{d}.variantQC.png",
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/donor{{d}}.labels.png",
               nclones=nclonelist, variants=params_clones["variants"]),
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/donor{{d}}.variants.labels.png",
               nclones=nclonelist, variants=params_clones["variants"])
    output:
        report("{out_dir}/results/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/donor{d}.pdf", category="Donor-specific")
    shell:
        "convert {input} {output[0]}"


rule mt_to_clones:
    input:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
                d=np.arange(config["N_DONORS"]))
    output:
        note="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/mtClones.ipynb",
        mtClones="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/cells_meta.tsv"
    params:
        don_dir = lambda wildcards, input: dirname(dirname(input[0])),
        outdir = lambda wildcards, output: dirname(output.note),
        rscript = join(ROOT_DIR, "src", "clones", "mt_to_clones.ipynb"),
        n_donors = config["N_DONORS"],
        samples = ",".join(samples['sample_name'].values)
    shell: "papermill -p don_dir {params.don_dir} -p outdir {params.outdir} -p t {wildcards.t} -p rt {wildcards.rt} -p n_donors {params.n_donors} -p samples {params.samples} {params.rscript} {output.note}"


########################################################################
## 5. Calling clones for each donor
########################################################################
rule vireo:
    input:
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb" #get_input_multiplex
    output:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{{nclones}}/donor{d}.variants.labels.png", d=np.arange(config["N_DONORS"])),#, category="lineage"),
        "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_simple/vireo/nclones{nclones}/cells_meta.tsv"
    params:
        notebook=join("src", "vireo", "2_MT_Lineage_Construct.ipynb"),
        #output_notebook = lambda wildcards, output: join(dirname(output[0]), 'clones.ipynb'),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["N_DONORS"],
        nclones= lambda wildcards: wildcards.nclones #",".join([str(x) for x in nclonelist])
    threads: 8
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p N_DONORS {params.N_DONORS} -p nclones {params.nclones} {params.notebook} {params.OUTDIR}/output.ipynb"


# ########################################################################
# ## Variants A: After multiplexing, separate by donor, grab variants from filters
# ## that overlap with both conditions, and then run mgatk to call variants again.
# ##
# # Extract pileups for each donor from original filter matrices
# ########################################################################
rule scPileup_filter_mgatk:
    input:
        cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
            sample=samples.index),
        mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb",
    output:
        cov = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{d}/d{d}.coverage.txt",
            d=range(config["N_DONORS"]))
    params:
        cells_meta = lambda wildcards, input: join(dirname(input.mult), "cells_meta.tsv"),
        sample = samples['sample_name'].values,
        #outdir = lambda wildcards, output: dirname(output.cov),
    script: join(ROOT_DIR, "src/donor_filter_mgatk.py")


module procDonMgatkMod:
    snakefile: "./rules/clone_detection/proc_donormgatk.smk"
    config: config

use rule * from procDonMgatkMod as donMGATK_*


# ########################################################################
# ## Variants B: After multiplexing, separate by donor, keep variants from filters
# ## that overlap with both conditions, and then use those for clonal calling (may include a lot of variants).
# ##
# # Extract pileups for each donor from original filter matrices
# ########################################################################
rule scPileup_init:
    input:
        cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
            sample=samples.index),
        mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        af = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_init/donor{d}/af.tsv",
                    d=np.arange(config["N_DONORS"]))
    params:
        sample = samples['sample_name'].values,
        cells_meta = lambda wildcards, input: join(dirname(input.mult),"cells_meta.tsv"),
        ref_mt = params_mgatk["chrM_refAllele"]
    script: join(ROOT_DIR, "src/donor_to_clonesInput/donor_to_af.py")


module procInitMod:
    snakefile: "./rules/clone_detection/proc_init.smk"
    config: config

use rule * from procInitMod as procInitMod_*


########################################################################
## Variants C: 'Simple': Directly from multiplex output.
## Variants not seen in one condition might be 0'd out, which can bias
## the clone calling. Only in Vireo
########################################################################

########################################################################
## Variants D: 'simpleUnion': Union of variants from both mgatk outputs. Similar to init,
# but it also adds in coverage if a variant was filtered out in a specific condition.
########################################################################
rule scPileup_simpleUnion:
    input:
        cov =  expand("{{out_dir}}/data/{sample}/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/{sample}.coverage.txt",
            sample=samples.index),
        mult = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    output:
        af = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/multiplex/clones_simpleUnion/donor{d}/af.tsv",
            d=np.arange(config["N_DONORS"]))
    params:
        sample = samples['sample_name'].values,
        cells_meta = lambda wildcards, input: join(dirname(input.mult),"cells_meta.tsv"),
        ref_mt = params_mgatk["chrM_refAllele"]
    script: join(ROOT_DIR, "src/donor_to_clonesInput/donor_to_af_simpleUnion.py")


module procSimUnMod:
    snakefile: "./rules/clone_detection/proc_simpleUnion.smk"
    config: config


use rule * from procSimUnMod as procSimUnMod_*

#########################################################################
## Clones stats and QC
#########################################################################
def get_counts_in(wildcards):
    w = wildcards
    print('w', w)
    if w.variants == "mgatkdonor":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/donor{w.d}/mgatk/d{w.d}.variant.rds"
    elif w.variants == "simple":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/multiplex.ipynb"
    elif w.variants == "init":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/clones_init/donor{w.d}/af.tsv"
    elif w.variants == "simpleUnion":
        return f"{w.out_dir}/data/merged/MT/cellr_{w.cellrbc}/numread_{w.num_read}/filters/minC{w.mincells}_minR{w.minreads}_topN{w.topN}_hetT{w.hetthresh}_hetC{w.minhetcells}_hetCount{w.hetcountthresh}_bq{w.bqthresh}/mgatk/vireoIn/multiplex/clones_simpleUnion/donor{w.d}/af.tsv"


module cloneMod:
    snakefile: "./rules/clone_detection/clone.smk"
    config: config

use rule * from cloneMod as clone_*


use rule convert_to_af from cloneMod as clone_convert_to_af with:
    input:
        counts_in = get_counts_in,
        cells_meta = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/sc_af/donor{d}/sc_af.ipynb",
    params:
        var_type = lambda wildcards: wildcards.variants,
        note = join(ROOT_DIR, "workflow", "notebooks", "clone_af_dendrograms", "Convert_to_AF.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note),
        indir = lambda wildcards, input: dirname(input.cells_meta),
        counts_indir = lambda wildcards, input: dirname(input.counts_in),


########################################################################
## Enrichment
########################################################################
module cloneEnrichMod:
    snakefile: "./rules/clone_detection/enrichment.smk"
    config: params

use rule * from cloneEnrichMod as cloneEnrich_*

############################################
## Compare methods:
############################################
rule compare_methods:
    input:
        expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
                variants=[x for x in params_clones["variants"] if x != "simple"], kparam=params_clones["knn"]["params"]["resolution"]),
        # expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv",
        #         variants=params_clones["variants"], nclones=nclonelist)
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/comparisons.ipynb",
        fig = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/",
                              "methods_nAgreeNorm.png", "methods_nAgreeNorm_agg.png", "methods_n_T_T_agg.png", "methods_nTTNorm_agg.png"), category="Methods"),
    params:
        all = lambda wildcards, input: ",".join(x for x in input),
        note = join("src", "clones_compare", "distance_matrix.ipynb"),
        outdir = lambda wildcards, output: dirname(output[0]),
        prefix = config["prefix"],
        cfg_outdir = config["outdir"],
        params_f = config["config"]
    shell:
        "papermill -p all_files {params.all} -p outdir {params.outdir} -p cfg_outdir {params.cfg_outdir} -p prefix {params.prefix} -p params_f {params.params_f} {params.note} {output.note}"

