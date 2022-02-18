from os.path import join, dirname
from src.config import ROOT_DIR

########################################################################################################################
# DE and GSEA rules
########################################################################################################################
def get_de_script(wildcards):
    if wildcards.de_group == "btwnChange_inClust":
        return join(ROOT_DIR,"workflow/notebooks/de_clones_change/DE_genes_btwnChange_inClust.ipynb")
    elif wildcards.de_group == "btwnChange":
        return join(ROOT_DIR,"workflow/notebooks/de_clones_change/DE_genes_btwnChange.ipynb")
    return

rule runDE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "SE.rds", #config["se_f"], #{outdir}/annotation_clones/SE.rds",
    output:
        note =  "{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/de.ipynb",
    params:
        rscript = get_de_script, #"workflow/notebooks/de_clones_change/DE_genes_btwnConds.ipynb",
        outdir = lambda wildcards, output: dirname(output.note),
        minPct=lambda wildcards: wildcards.minPct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        samples = ",".join(config["samples"].index),
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p sample_names {params.samples} {params.rscript} {output.note}"


rule summarizeDE:
    input:
        note =  "{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/output_de.ipynb",
    output:
        note =  "{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/output_summarizeDE_{pthresh}.ipynb",
    params:
        indir = lambda wildcards, input: dirname(input.note),
        outdir = lambda wildcards, output: dirname(output.note),
        rscript = get_de_script, #"workflow/notebooks/de_clones_change/DE_genes_btwnConds.ipynb",
        top_de=3,
        pthresh=lambda wildcards: wildcards.pthresh,
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p top_de {params.top_de} -p pthresh {wildcards.pthresh} {params.rscript} {output.note}"


rule runGSEA:
    input:
        DE_out_path ="{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/de.ipynb"
    output:
        note= "{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/GSEA_pthresh.{gsea_pval}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input.DE_out_path)),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/runGSEA_btwnCond.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule runGSEA_compact:
    input:
        DE_out_path ="{de_indir}/de.ipynb"
    output:
        note= "{de_indir}/GSEA_pthresh.{gsea_pval}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input.DE_out_path)),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/runGSEA_btwnCond.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule runGSEA_fullParams:
    input:
        DE_out_path ="{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/de.ipynb"
    output:
        note= "{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_p.{pthresh}__stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(f"{dirname(input[0])}_pthresh_{wildcards.pthresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/runGSEA_btwnCond.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval,
        stat_col = lambda wildcards: wildcards.stat_col,
        prefilter = False, #lambda wildcards: wildcards.prefilter,
        padjmethod = lambda wildcards: wildcards.padjmethod,
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p padjmethod {params.padjmethod} -p pthresh {wildcards.pthresh} -p gsea_pthresh {params.gsea_pval} -p prefilter {params.prefilter} -p stat_col {params.stat_col:q} -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"


rule summaryGSEA:
    input:
        note= "{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/GSEA_pthresh.{gsea_pval}/GSEA.ipynb",
    output:
        note= "{outdir}/de/{de_group}/minPct_{minPct}_logfc{logfc_threshold}/GSEA_pthresh.{gsea_pval}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input.note), #join(f"{dirname(input[0])}_pthresh_{wildcards.pthresh}", "btwnConds_inClust"),
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"


rule summaryGSEA_compact:
    input:
        "{gsea_dir}/GSEA.ipynb",
    output:
        note="{gsea_dir}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.pthresh}", "btwnConds_inClust"),
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"

