from os.path import join, dirname
from src.config import ROOT_DIR

########################################################################################################################
# Run DE and GSEA in expanded vs regressed clones for genes and peaks
# Requires clones file with "change" column as either "expansion", "regression", "no_change"
########################################################################################################################
rule preprocess:
    input:
        enrichment_dir = "enrichmentNorm.csv"
    output:
        clone_change_f="{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
        note="{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/preproc.ipynb"
    params:
        outdir = lambda wildcards, output: dirname(output.note),
        script = "workflow/notebooks/de_clones_change/preprocess_clones.ipynb"
    shell:
        "papermill -p clone_change_f {output.clone_change_f} -p clone_filt {wildcards.filt} -p outdir {params.outdir} "
        "-p min_cell {wildcards.min_cell} -p min_cell_both {wildcards.min_cell_both} "
        "-p pthresh {wildcards.pthresh} -p shuffle_stats {wildcards.shuffle} "
        "{params.script} {output.note} "



rule btwnClone_DE:
    """ Compares clusters to each other. For now works with Gene Activity
    
    Not tested with peaks as gene activity yet.
    """
    input:
        se_f = "SE.rds", #config["se_f"], #{outdir}/annotation_clones/SE.rds",
        clone_change_f = "{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/preproc/clone_change.csv",
    output:
        note =  "{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/de.ipynb",
    params:
        rscript = "workflow/notebooks/de_clones_change/DE_genes_btwnConds.ipynb",
        outdir = lambda wildcards, output: dirname(output.note),
        minPct=lambda wildcards: wildcards.minPct,
        logfcthresh= lambda wildcards: wildcards.logfc_threshold,
        top_de=3,
        samples = ",".join(config["samples"].index),
        p_thresh=lambda wildcards: wildcards.p_thresh,
    shell: "papermill -p se_f {input.se_f} -p clone_change_f {input.clone_change_f} -p outdir {params.outdir} -p top_de {params.top_de} -p sample_names {params.samples} {params.rscript} {output.note}"


rule btwnClone_runGSEA_sepDonors:
    input:
        note = "{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/de.ipynb",
    output:
        note= "{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "sepDonors"), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = "workflow/notebooks/de_clones_change/runGSEA.ipynb",
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        #pattern = "donor_*",
        #is_dir="FALSE"
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p padjmethod {wildcards.padjmethod} -p pthresh {wildcards.p_thresh} -p gsea_pthresh {wildcards.gsea_pval} -p prefilter {wildcards.prefilter} -p stat_col {wildcards.stat_col:q} -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"


rule btwnClone_summaryGSEA_sepDonors:
    input:
        note="{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb"
    output:
        note="{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/sepDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]),
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"


rule btwnClone_runGSEA_allDonors:
    input:
        note = "{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/de.ipynb",
    output:
        note= "{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/allDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "sepDonors"), #join(f"{dirname(input[0])}_pthresh_{wildcards.p_thresh}", "btwnConds_inClust"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = "workflow/notebooks/de_clones_change/runGSEA.ipynb",
        gsea_dir = join(ROOT_DIR, "software/Bioinformatics_Tools/"), #/ is necessary
        #pattern = "donor_*",
        #is_dir="FALSE"
    #conda_env = "../envs/gsea_manual.yml"
    conda:
        "../envs/gsea_manual.yml" #environment with clusterprofiler
    shell: "papermill -p padjmethod {wildcards.padjmethod} -p pthresh {wildcards.p_thresh} -p gsea_pthresh {wildcards.gsea_pval} -p prefilter {wildcards.prefilter} -p stat_col {wildcards.stat_col:q} -p DE.out.path {params.input} -p export.path {params.output} -p gsea_dir {params.gsea_dir} {params.rscript} {output[0]}"


rule btwnClone_summaryGSEA_allDonors:
    input:
        "{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/allDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/GSEA.ipynb"
    output:
        note="{outdir}/clones_change/filt_{filt}__shuffle_{shuffle}__padj_{use_p_adjust}__pthresh_{pthresh}_minC_{min_cell}__bothMinC__{min_cell_both}/de/minPct_{minPct}_logfc{logfc_threshold}_pthresh_{p_thresh}/allDonors/GSEA_pthresh.{gsea_pval}_pref.{prefilter}_stat.{stat_col}_padj.{padjmethod}/summary.ipynb",
    params:
        input = lambda wildcards, input: dirname(input[0]),
        rscript = join(ROOT_DIR, "R_scripts/annotation_clones/summarizeGSEA.ipynb"),
    shell: "papermill -p export_path {params.input} {params.rscript} {output.note}"

