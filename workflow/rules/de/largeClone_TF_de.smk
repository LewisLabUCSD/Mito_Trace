from os.path import join, dirname
from src.config import ROOT_DIR
########################################################################################################################
## TF analysis on larger clones
########################################################################################################################
rule setup_motifs_largeClones:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
    output:
        note="{outdir}/annotation_clones/DE_large/output.ipynb",
        largeclones_se = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    threads: 16
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        n_donors = config["N_DONORS"],
        rscript= join(ROOT_DIR, "R_scripts/annotation_clones/de_TF_largeClones/setup_motifs_large_clones.ipynb")
    shell: "papermill -p se_f {input.se_f} -p outdir {params.outdir} -p n_donors {params.n_donors} {params.rscript} {output.note}"


rule runDE_TF_largeClones:
    input:
        se_f = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        note = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        out2 = report(directory(expand("{{outdir}}/annotation_clones/DE_large/minPct_{{minPct}}__logThresh_{{logThresh}}/donor{d}_TF",
            d = range(config["N_DONORS"]))),
            patterns=["*png", "*csv"],
            category="lineage", subcategory="donor {d}"),
    threads: 8
    params:
        indir = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        n_donors = config["N_DONORS"],
        rscript= join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_TF_largeClones/DE_clones.TF.vCurrent.ipynb"),
    #min_pct = lambda wildcards
    shell: "papermill -p indir {params.indir} -p outdir {params.outdir} -p min_pct {wildcards.minPct} -p logfc_thresh {wildcards.logThresh} -p n_donors {params.n_donors} {params.rscript} {output.note}"


rule summary_TF_largeClones:
    input:
        de = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        out = multiext("{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
            "output_Summary.ipynb", "dotplot.allDonors.top.pdf", "dotplot.allDonors.clones.top.pdf",
            "heatmap.allDonors.top.pdf",
            "heatmap.allDonors.split.top.pdf", "large.clones.cdf.png"),
        donorDot = expand("{{outdir}}/annotation_clones/DE_large/minPct_{{minPct}}__logThresh_{{logThresh}}/cdf_thresh__{{cdf_thresh}}/donor{d}_TF/dot.top.png",
            d=range(config["N_DONORS"])),
    #note="{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/output_Summary.ipynb",
    params:
        indir = lambda wildcards, input: dirname(input.de),
        se_indir = lambda wildcards, input: dirname(input.se_f),
        n_donors = config["N_DONORS"],
        cdf_thresh = lambda wildcards: wildcards.cdf_thresh,
        rscript= join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_TF_largeClones/DE_clones.TF.summaryPlot.ipynb")
    shell: "papermill -p indir {params.indir} -p se_indir {params.se_indir} -p n_donors {params.n_donors} -p cdf_thresh {params.cdf_thresh} {params.rscript} {output[0]}" #{output[0]}

rule output_largeclones:
    """ Summarizes TF results across all clonal comparisons. 
    """
    input:
        de = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/output_DE.ipynb",
        se_f = "{outdir}/annotation_clones/DE_large/se.clonesfilt.rds"
    output:
        note = "{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/TFactivity.ipynb",
        out = multiext("{outdir}/annotation_clones/DE_large/minPct_{minPct}__logThresh_{logThresh}/cdf_thresh__{cdf_thresh}/",
            "dotplot.allDonors.clones.topPvals.pdf", "dotplot.allDonors.clones.topPvals.noZ.pdf",
            "clone.TFactivity.topPvals.txt")
    params:
        indir = lambda wildcards, input: dirname(input.de),
        se_indir = lambda wildcards, input: dirname(input.se_f),
        n_donors = config["N_DONORS"],
        cdf_thresh = lambda wildcards: wildcards.cdf_thresh,
        rscript = join(ROOT_DIR,"workflow/notebooks/nuclear_de/de_TF_largeClones/DE_clones.TF.summaryPlot.DonorConcat.Development.ipynb")
    shell: "papermill -p indir {params.indir} -p se_indir {params.se_indir} -p n_donors {params.n_donors} -p cdf_thresh {params.cdf_thresh} {params.rscript} {output.note}"

rule runGSEA:
    """ GSEA on large clones DE
            
        TODO: implement
    """
    input:
        DE_out_path = "{outdir}/annotation_clones/DE/donor{n}/DE.ipynb"
    output:
        "{outdir}/annotation_clones/DE/donor{n}/GSEA/clusters/GSEA.ipynb",
    params:
        input = lambda wildcards, input: join(dirname(input[0]), "clusters"),
        output = lambda wildcards, output: dirname(output[0])+"/", #/ needed for the GSEA script
        rscript = join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_TF_largeClones/clones/runGSEA_clones.ipynb")
    shell: "papermill -p DE.out.path {params.input} -p export.path {params.output} {params.rscript} {output[0]}"



rule runDE_enrichedClones:
    input:
        se_f = "{outdir}/annotation_clones/SE.rds",
        enrich_f = expand("{{outdir}}/enrichmentNorm_donor{d}.csv",
            d=range(config["N_DONORS"])),
    output:
        note="{outdir}/annotation_clones/DE/DE.ipynb"
    params:
        #d = lambda wildcards: wildcards.d,
        enrich_d = lambda wildcards, input: dirname(input[0]),
        outdir = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["N_DONORS"],
        rscript= join(ROOT_DIR, "workflow/notebooks/nuclear_de/de_TF_largeClones/addEnrichment.vCurrent.ipynb")
    shell: "papermill -p enrich_d {params.enrich_d} -p se_f {input.se_f} -p outdir {params.outdir} -p N_DONORS {params.N_DONORS} {params.rscript} {output.note}"
