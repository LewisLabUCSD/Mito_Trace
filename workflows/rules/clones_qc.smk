from os.path import join, dirname
################################
## Clones QC: Transitions+#distinguishing variants per clone
################################
rule distinguishing_variants:
    input: "{outdir}/cells_meta.tsv"
    params:
        config["clones_qc"]["params"]["binary_thresh"]
    output: expand("{{outdir}}/clones_qc/{{params}}/{out}",
                   out=["transitions.png",
                        "distinguishing_variants.png"])



################################
## MT Heteroplamsy across conditions
################################
rule condition_heteroplasmy:
    input: "{outdir}/cells_meta.tsv"
    output:
        note = "{outdir}/clones_qc/{params}/output.ipynb",
        het = "{outdir}/clones_qc/{params}/heteroplasmy.png",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        note = join("src", "clones", "condition_heteroplasmy.ipynb")
    shell: "papermill -p INDIR {params.INDIR} {params.note} {output}"


rule condition_shuffle_heteroplasmy:
    input: "cells_meta.tsv"
    output:
        note = "{outdir}/clones_qc/{params}/shuffle_het/output.ipynb",
        dist = "{outdir}/clones_qc/{params}/shuffle_het/distribution.png",
        vars_sig = "{outdir}/clones_qc/{params}/shuffle_het/vars_sig.txt",
        vars_sig_het = "{outdir}/clones_qc/{params}/shuffle_het/vars_sig_het.png",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        note = join("src", "clones", "shuffle_het.ipynb")
    shell: "papermill -p INDIR {params.INDIR} {params.note} {output}"


################################
## Enrichment QC
################################
rule shuffle_enrichment:
    input:
        clones_f = "{outdir}/cells_meta.tsv",
        enrich_f = "{outdir}/enrichment/enrichment.csv"
    output:
        note = "{outdir}/clones_qc/{params}/shuffle_enrich/output.ipynb",
        dist = "{outdir}/clones_qc/{params}/shuffle_enrich/distribution.png",
        vars_sig = "{outdir}/clones_qc/{params}/shuffle_enrich/clones_sig.txt",
        #vars_sig_het = "{outdir}/clones_qc/{params}/shuffle_enrich/clones_sig_het.png",
    params:
        INDIR = lambda wildcards, input: dirname(input[0]),
        note = join("src", "clones", "shuffle_enrich.ipynb")
    shell: "papermill -p INDIR {params.INDIR} {params.note} {output}"



