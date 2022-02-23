from os.path import join, dirname
from src.config import ROOT_DIR


## inClone_btwnCond
rule inClone_btwnCond_DE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note =  "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    params:
        rscript= join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_conditions/DE_genes_inClones_btwnConditions.ipynb"), # The script defaults to the granja data
        samples = ",".join(config["samples"].index),
        outdir = lambda wildcards, output: dirname(output.note), #"{outdir}/annotation_clones/de_btwnConds_{assay}/minPct_{btwnMinpct}_logfc{logfcthresh}/pthresh_{p_thresh}",
        minPct=lambda wildcards: wildcards.btwnMinpct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        p_thresh=lambda wildcards: wildcards.p_thresh,
    # test_use="wilcox",
    # latent_vars="NULL",
    threads: 8
    shell: "papermill -p se_f {input.se_f} -p min_cells {wildcards.min_cells} -p p_thresh {params.p_thresh} -p outdir {params.outdir} -p top_de {params.top_de} -p sample_names {params.samples} {params.rscript} {output.note}"

rule runGSEA_inClone_btwnCond:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_conditions/runGSEA_inClones_btwnConditions.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval,
        stat_col = lambda wildcards: wildcards.stat_col,
        prefilter = lambda wildcards: wildcards.prefilter,
        padjmethod = lambda wildcards: wildcards.padjmethod,
        pthresh = lambda wildcards: wildcards.p_thresh,
        pattern = "*clone_.*counts.*"
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill  -p DE.out.path {params.input} -p export.path {params.output} -p pattern {params.pattern:q} -p padjmethod {params.padjmethod} -p pthresh {params.pthresh} -p gsea_pthresh {params.gsea_pval} -p prefilter {params.prefilter} -p stat_col {params.stat_col:q} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule summaryGSEA_inClone_btwnCond:
    input:
        "{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    output:
        note="{outdir}/annotation_clones/de_clone_btwncond_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        script = join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_conditions/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.script} {output.note}"

####################
## Between Variants!
####################
rule inVar_btwnCond_DE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
        mt_cells = "cells_meta.tsv" if "cells_meta" not in config else config["cells_meta"]
    output:
        note =  "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    params:
        rscript= join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_conditions/DE_genes_inVariants_btwnConditions.ipynb"), # The script defaults to the granja data
        samples = ",".join(config["samples"].index),
        outdir = lambda wildcards, output: dirname(output.note), #f"{dirname(output.note)}_pthresh_{wildcards.p_thresh}", #"{outdir}/annotation_clones/de_btwnConds_{assay}/minPct_{btwnMinpct}_logfc{logfcthresh}/pthresh_{p_thresh}",
        minPct=lambda wildcards: wildcards.btwnMinpct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        p_thresh=lambda wildcards: wildcards.p_thresh,
    # test_use="wilcox",
    # latent_vars="NULL",
    threads: 8
    shell: "papermill -p se_f {input.se_f} -p cells_meta_f {input.mt_cells} -p outdir {params.outdir} -p min_cells {wildcards.min_cells} -p p_thresh {params.p_thresh} -p top_de {params.top_de} -p sample_names {params.samples} {params.rscript} {output.note}"

rule runGSEA_inVars_btwnCond:
    input:
        DE_out_path = "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/de.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_conditions/runGSEA_inClones_btwnConditions.ipynb"),
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        gsea_pval = lambda wildcards: wildcards.gsea_pval,
        stat_col = lambda wildcards: wildcards.stat_col,
        prefilter = lambda wildcards: wildcards.prefilter,
        padjmethod = lambda wildcards: wildcards.padjmethod,
        pthresh = lambda wildcards: wildcards.p_thresh,
        pattern = "*clone_*"
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p pattern {params.pattern:q} -p padjmethod {params.padjmethod} -p pthresh {params.pthresh} -p gsea_pthresh {params.gsea_pval} -p prefilter {params.prefilter} -p stat_col {params.stat_col:q} -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"

rule summaryGSEA_inVars_btwnCond:
    input:
        "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    output:
        note= "{outdir}/annotation_clones/de_clone_btwnvars_RNA/minPct_{btwnMinpct}_logfc{logfc_threshold}_minCells{min_cells}_pthresh{p_thresh}/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        rscript = join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_conditions/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"

