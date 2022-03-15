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


########################################################################
# Setup parameters and outdir
########################################################################
res = join(config["outdir"], "pipeline", config["prefix"])

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


params["markers_f"] = config["markers_f"]

# Add cell markers to annotation:
cluster_names_f = join(res, "results", "atac_clusters.txt")
if not os.path.exists(cluster_names_f):
    cluster_names_f = "None"
params["cluster_names_f"] = cluster_names_f
#workdir: config["work_dir"]


###################################
## Rules and workflow
###################################
###################################
rule all:
    input:
        # A. Preproc: Initial MT coverage per position per cell
        expand("{out_dir}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png",
            out_dir=res,sample=samples.index, num_read=num_reads_filter, cellrbc=cellrbc),

        # B. Multiplex: Multiplexed VAF and depth dendrograms for each donor and selected variants
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

        # Annotation clones
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/_nuclear_clones_complete.txt",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            variants=[x for x in params_clones["variants"] if x != "simple"],
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'], gff=gff,
            kparam=params_clones["knn"]["params"]["resolution"]),

        # MT Plots
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA_af/mtPlots.ipynb",
            out_dir=res,cellrbc=cellrbc,num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'],bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            variants=[x for x in params_clones["variants"] if x != "simple"], gff=gff,
            t=params_annclo["t"], rt=params_annclo["rt"]),


        ##  Methods compare
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/comparisons/comparisons.ipynb",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh']),

        # Merge clone enrich and nuclear
        expand("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage/out.ipynb",
            out_dir=res,cellrbc=cellrbc, num_read=num_reads_filter,
            mincells=ft['mincells'],minreads=ft['minreads'],topN=ft["topN"],hetthresh=ft['hetthresh'],minhetcells=ft['minhetcells'],
            hetcountthresh=ft['hetcountthresh'], bqthresh=ft['bqthresh'],
            kparam=config["params_clones"]["knn"]["params"]["resolution"],
            minPct=params_annclo["min_pct"], gsea_pval=params_annclo["gsea_pval"], #assay=params_annclo["assay"],
            gff=gff, min_cells=params_annclo["min_cells"],
            variants=[x for x in params_clones["variants"] if x != "simple"],
            stat_col=params_annclo["stat_col"], padjmethod=params_annclo["padjmethod"],
            filt = params_clch["filt"], shuffle=params_clch["shuffle"], use_p_adjust=params_clch["use_p_adjust"],
            pthresh=params_clch["pthresh"], min_cell=params_clch["min_cell"], min_cell_both=params_clch["min_cell_both"],),

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
########################################################################


########################################################################
## 1. Go from 10x output to filtered scPileup outdir.
########################################################################
use rule * from mtpreprocMod

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


# def check_conds(wildcards):
#     w = wildcards
#     if samples.shape[0] == 1:
#         return f"{w.out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/cellSNP.tag.AD.mtx"
#
# rule one_cond_merged:
#     input:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/cellSNP.tag.AD.mtx"
#     output:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/cellSNP.tag.AD.mtx"

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

# use rule knn_mgatkdonor_concat from procDonMgatkMod as donMGATK_knn_mgatkdonor_concat with:
#     input:
#         cells = expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{{kparam}}/labels/donor{d}/cells_meta.tsv",
#                 d = np.arange(config["N_DONORS"])),
#     output:
#         cells_all = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/cells_meta.tsv",
# #         expand("{{outdir}}/clones/variants_mgatkdonor/vireo/nclones{{nclones}}/donor{d}_cells_meta.tsv",
# #                d = np.arange(config["N_DONORS"])),

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


# use rule barcodes_inClones from cloneMod as clone_barcodes_inClones with:
#     input:
#         counts_in = get_counts_in,
#         cells_meta = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/barcodes/btwnClones/donor{d}.ipynb",
#     params:
#         var_type = lambda wildcards: wildcards.variants,
#         note=join(ROOT_DIR, "worklow/notebooks/clone_af_dendrograms", "MT_barcodes/inClones.ipynb"),
#         outdir=lambda wildcards, output: dirname(output.note),
#         indir= lambda wildcards, input: dirname(input.cells_meta),
#         counts_indir = lambda wildcards, input: dirname(input.counts_in)


########################################################################
## Enrichment
########################################################################
module cloneEnrichMod:
    snakefile: "./rules/clone_detection/enrichment.smk"
    config: params

use rule * from cloneEnrichMod as cloneEnrich_*

# rule vireo_enrichment:
#     version: '1.0' # no +1 norm error
#     input:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv"
#     output:
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png",
#             category="lineage", subcategory="enrichment"),
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/enrichmentNorm.csv",
#             category="lineage", subcategory="enrichment")
#     params:
#         #clones_indir = lambda wildcards, input: dirname(input[0]),#lambda wildcards, input: dirname(input[0]),
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         script = join("src", "lineage", "lineage_enrichment.py"),
#         samples=",".join(samples.index),
#         comparisons = config["comparisons"] if "comparisons" in config else "None"
#     shell: "python {params.script} {input} {params.OUTDIR} {params.samples} --tests {params.comparisons:q}"
#
#
#
# rule knn_enrichment:
#     input:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv",
#     output:
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png",
#             category="lineage", subcategory="enrichment"),
#         report("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/enrichmentNorm.csv",
#             category="lineage", subcategory="enrichment")
#     params:
#         OUTDIR = lambda wildcards, output: dirname(output[0]),
#         script = join("src", "lineage", "lineage_enrichment.py"),
#         samples = ",".join(samples.index),
#         comparisons = config["comparisons"] if "comparisons" in config else "None"
#     shell: "python {params.script} {input} {params.OUTDIR} {params.samples} --tests {params.comparisons:q}"




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


############################################
## Annotation
############################################
module annoMod:
    snakefile: "./rules/annotation_pipe.smk"
    config: config

use rule * from annoMod as anno_*

use rule integrateSignac from annoMod as anno_integrateSignac with:
    input:
       a = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.rds",
    output:
        a = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.ipynb",
        b = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        c = report(expand("{{out_dir}}/data/annotation/gff_{{gff}}/mergedSamples/{f}",
                      f=["QC_02.png",
                         "integrated.merged.compare.png",
                         "integrated.batch.png",
                         "integrated.lsi.clusters.png"]), category="Nuclear ATAC", subcategory="UMAP")


############################################
## Clones DE
############################################
module annCloMod:
    snakefile: "./rules/annotation_clones.smk"
    config: params

use rule * from annCloMod as annClo_*

use rule addClones from annCloMod as annClo_addClones with:
    input:
        noc = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        clones = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/cells_meta.tsv"
    output:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/SE.rds",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/addClones.ipynb"

use rule addClones from annCloMod as knn_annClo_addClones with:
    input:
        noc = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds",
        clones = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/cells_meta.tsv"
    output:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/addClones.ipynb"


############################################
## clone expansion vs regressed DE
############################################
module cloneChangeMod:
    snakefile: "./rules/de/de_clones_change.smk"
    config: params #params_clch

use rule * from cloneChangeMod as cloneChange_*

use rule preprocess from cloneChangeMod as cloneChange_preprocess with:
    input:
        enrichment_dir = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/enrichmentNorm.csv"
    output:
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/preproc.ipynb"


## get clone meta from stats and immune clusters
rule merge_enrich_and_lineage:
    input:
        clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        dominant_cluster = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/combinedConditions_dominant_cluster.csv"
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage/out.ipynb",
    params:
        dom_indir = lambda wildcards, input: dirname(input.dominant_cluster),
        script = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_count_embed", "aggregate_clone_meta.ipynb"),
        outdir = lambda wildcards, output: dirname(output.note)
    shell:
        "papermill -p dom_indir {params.dom_indir} -p enrichment_f {input.clone_change_f} -p outdir {params.outdir} {params.script} {output.note}"

# rule merge_enrich_and_lineage_and_input:
#     input:
#         clone_change_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
#         dominant_cluster = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/combinedConditions_dominant_cluster.csv",
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/merge_enrich_and_lineage_and_input/out.ipynb",
#     params:
#         dom_indir = lambda wildcards, input: dirname(input.dominant_cluster),
#         script = join(ROOT_DIR, "workflow/notebooks/lineage_clones/clone_cluster_count_embed", "aggregate_clone_meta.ipynb"),
#         outdir = lambda wildcards, output: dirname(output.note)
#     shell:
#         "papermill -p se_f {input.se_f} -p dom_indir {params.dom_indir} -p enrichment_f {input.clone_change_f} -p outdir {params.outdir} {params.script} {output.note}"


# Overlay mt variants on embeddings
rule mtVarsPlot:
    input:
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
        mt_cells = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/cells_meta.tsv"
    output:
        note =  "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA_af/mtPlots.ipynb",
    params:
        mt_cells = lambda wildcards, input: input.mt_cells.replace("cells_meta.tsv", "cells_mt.tsv"),
        rscript= join(ROOT_DIR, "workflow/notebooks/lineage_clones//mtVarsPlot.ipynb"), # The script defaults to the granja data
        outdir = lambda wildcards, output: dirname(output.note),
    # test_use="wilcox",
    # latent_vars="NULL",
    threads: 8
    shell: "papermill -p se_f {input.se_f} -p cells_meta_f {input.mt_cells} -p outdir {params.outdir} {params.rscript} {output.note}"

## DE without clone
module btwnClustersMod:
    snakefile: "./rules/de/btwnClusters.smk"
    config: params
use rule * from btwnClustersMod as btwnClusters_*

module btwnConditionsMod:
    snakefile: "./rules/de/btwnConditions.smk"
    config: params
use rule * from btwnConditionsMod as btwnConditionsMod_*

## DE with clone
module inCloneBtwnConditionsMod:
    snakefile: "./rules/de/inClones_btwnConditions.smk"
    config: params

use rule * from inCloneBtwnConditionsMod as inCloneBtwnConditionsMod_*



use rule inVar_btwnCond_DE from annCloMod as annClo_inVar_btwnCond_DE with:
    input:
        mt_cells = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/multiplex/mt_clones/thr__{t}_rt__{rt}/cells_meta.tsv",
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/SE.rds",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/thr__{t}_rt__{rt}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",



# use rule dominant_clone_clust_in_clone from annCloMod as annClo_dominant_clone_clust_in_clone with:
#     params:
#         cluster_names_f = cluster_names_f
    # input:
    #     se_meta="{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/se_cells_meta.tsv",
    #     cluster_names_f = lambda wildcards: get_clust_markers
    # output:
    #     note= "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/dominant.ipynb",
    #     results = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/dominant_clone_clust/",
    #         "dominant_cluster.csv", "cluster_clone_meta.csv", "cluster_clone_counts_normalized.csv",
    #         "combinedConditions_dominant_cluster.csv", "combinedConditions_cluster_clone_meta.csv", "combinedConditions_cluster_clone_counts_normalized.csv"))

## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb"
#

# ## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as knn_annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb",


## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         out = multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
#                        "output_Summary.ipynb", "dotplot.allDonors.top.pdf", "dotplot.allDonors.clones.top.pdf",
#                        "heatmap.allDonors.split.top.pdf", "large.clones.cdf.png"),
#
# ## Use annotation_clones to run DE on clones
# use rule summary_TF_largeClones from annCloMod as knn_annClo_summary_TF_largeClones with:
#     input:
#         se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
#         de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb"
#     output:
#         out = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
#                          "output_Summary.ipynb", "dotplot.allDonors.top.pdf", "dotplot.allDonors.clones.top.pdf",
#                          "heatmap.allDonors.split.top.pdf", "large.clones.cdf.png"),
#                          category="lineage and nuclear", subcategory="Method: knn Variants: {variants} resolution: {kparam} minPct: {minPct} logThresh: {logThresh}"),
#         donorDot = report(expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{{variants}}/knn/kparam_{{kparam}}/gff_{{gff}}/annotation_clones/DE_large/minPct_{{minPct}}__logThresh_{{logThresh}}/cdf_thresh__{{cdf_thresh}}/donor{d}_TF/dot.top.png",
#                                  d=range(config["N_DONORS"])),
#                           category="lineage and nuclear", subcategory="Method: KNN Variants: {variants} k resolution: {kparam} minPct: {minPct} logThresh: {logThresh}")


############################################
############################################
rule clones_donors_table:
    input:
        cells_meta_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/cells_meta.tsv"
    output:
        summary_table = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/summary/counts_summary.csv",
        summary_fig = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/nclones{nclones}/summary/counts_summary.pdf"
    params:
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        notebook = join(ROOT_DIR, "src", "clones")
    shell: "papermill -p cells_meta_f {input.cells_meta_f} -p OUTDIR {params.OUTDIR} {params.notebook} {params.OUTDIR}/output.ipynb"


use rule output_largeclones from annCloMod as annClo_output_largeclones with:
    input:
        de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds",
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/TFactivity.ipynb",
        out = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
            "dotplot.allDonors.clones.topPvals.pdf","dotplot.allDonors.clones.topPvals.noZ.pdf",
            "clone.TFactivity.topPvals.txt"),
            category="lineage and nuclear", subcategory="Method: vireo Variants: {variants} nClones: {nclones} minPct: {minPct} logThresh: {logThresh}")

## Use annotation_clones to run DE on clones
use rule output_largeclones from annCloMod as knn_annClo_output_largeclones with:
    input:
        de = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/TFactivity.ipynb",
        out = report(multiext("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
            "dotplot.allDonors.clones.topPvals.pdf","dotplot.allDonors.clones.topPvals.noZ.pdf",
            "clone.TFactivity.topPvals.txt"),
            category="lineage and nuclear", subcategory="Method: knn Variants: {variants} resolution: {kparam} minPct: {minPct} logThresh: {logThresh}"),




use rule runDE_TF from annoMod as anno_runDE_TF with:
    input:
        integrated="{out_dir}/data/annotation/gff_{gff}/mergedSamples/allSamples.integrated.rds"
    output:
        note="{out_dir}/data/annotation/gff_{gff}/mergedSamples/DE_TF/DE_TF.ipynb"


use rule runGSEA from annoMod as anno_runGSEA with:
    input:
        DE_out_path = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/DE/DE.ipynb"
    output:
        note = "{out_dir}/data/annotation/gff_{gff}/mergedSamples/DE/GSEA/clusters/GSEA.ipynb"

# use rule btwnClone_runGSEA_sepDonors from cloneChangeMod as cloneChange_btwnClone_runGSEA_sepDonors with:
#     input:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/de.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb"
#
#
# use rule btwnClone_summaryGSEA_sepDonors from cloneChangeMod as cloneChange_btwnClone_summaryGSEA_sepDonors with:
#     input:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb"
#     output:
#         note = "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
#         #note2= "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb"
#

# Summarize over clone results
# rule vireo_process:
#     input:
#         expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/volcano_Fisher_foldNorm.png",
#                 nclones=nclonelist,
#                 variants=params_clones["variants"]),
#         expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb",
#                 nclones=nclonelist,
#                 variants=params_clones["variants"],
#                 logThresh=params["annotation_clones"]["params"]["logfc_threshold"],
#                 minPct=params["annotation_clones"]["params"]["min_pct"],
#                 cdf_thresh=params["annotation_clones"]["params"]["cdf_thresh"],
#                 gff = [gff] ),
#         expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/enrichment/shuffle_stats/shuffle_stats.csv",
#                 nclones=nclonelist, variants=params_clones["variants"]),
#         expand("{{out_dir}}/results/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/vireo/nclones{nclones}/clones_dendro/donor{d}.{f}na.clust.max2.AF.png",
#                 d=np.arange(config["N_DONORS"]), nclones=nclonelist,
#                 variants=params_clones["variants"], f=["", "NoCondition."]),
#     output:
#         temp("{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/temp/vireo/.tmp")
#     shell: "touch {out_dir}"
#
#
# rule knn_process:
#     input:
#         expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png",
#                 kparam=config["params_clones"]["knn"]["params"]["resolution"],
#                 variants=[x for x in params_clones["variants"] if x != "simple"]),
#         # expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/gff_{gff}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb",
#         #         kparam=config["params_clones"]["knn"]["params"]["resolution"],
#         #         variants=[x for x in params_clones["variants"] if x != "simple"],
#         #         logThresh=params["annotation_clones"]["params"]["logfc_threshold"],
#         #         minPct=params["annotation_clones"]["params"]["min_pct"],
#         #         cdf_thresh=params["annotation_clones"]["params"]["cdf_thresh"],
#         #         gff=[gff]),
#         expand("{{out_dir}}/data/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/enrichment/shuffle_stats/shuffle_stats.csv",
#                 kparam=config["params_clones"]["knn"]["params"]["resolution"],
#                 variants=[x for x in params_clones["variants"] if x != "simple"]),
#         expand("{{out_dir}}/results/merged/MT/cellr_{{cellrbc}}/numread_{{num_read}}/filters/minC{{mincells}}_minR{{minreads}}_topN{{topN}}_hetT{{hetthresh}}_hetC{{minhetcells}}_hetCount{{hetcountthresh}}_bq{{bqthresh}}/mgatk/vireoIn/clones/variants_{variants}/knn/kparam_{kparam}/clones_dendro/donor{d}.{f}na.clust.max2.AF.png",
#                 d=np.arange(config["N_DONORS"]),kparam=config["params_clones"]["knn"]["params"]["resolution"],
#                 f=["", "NoCondition."], variants=[x for x in params_clones["variants"] if x != "simple"]),
#     output:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/temp/knn/.tmp"#.pipeline"
#     shell:
#         "touch {out_dir}"
#
#
#
# rule complete_lineage:
#     input:
#         "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/temp/{method}/.tmp",
#     output: "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/{method}/.completed"
#     shell: "touch {out_dir}"

#########################################################################
#########################################################################

# rule terminate:
#     input:
#          rules.complete_lineage.output
#     output:
#           "{out_dir}/data/merged/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_{variants}/{method}/.pipeline"
#     shell: "touch {out_dir}"
#
#
