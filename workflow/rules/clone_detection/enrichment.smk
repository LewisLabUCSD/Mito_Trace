from os.path import join, dirname

rule enrichment:
    input:
        "{outdir}/cells_meta.tsv",
    output:
        report("{outdir}/enrichment/volcano_Fisher_foldNorm.png",
            category="lineage", subcategory="enrichment"),
        report("{outdir}/enrichment/enrichmentNorm.csv",
            category="lineage", subcategory="enrichment")
    params:
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        script = join("src", "lineage", "lineage_enrichment.py"),
        samples = ",".join(config["samples"].index),
        comparisons = config["comparisons"] if "comparisons" in config else "None"
    shell: "python {params.script} {input} {params.OUTDIR} {params.samples} --tests {params.comparisons:q}"


rule clone_shuffle_stats:
    input:
        enrich_f = "{outdir}/enrichment/enrichmentNorm.csv",
        cells_meta_f = "{outdir}/cells_meta.tsv",
    output:
        "{outdir}/enrichment/shuffle_stats/output.ipynb",
        report("{outdir}/enrichment/shuffle_stats/shuffle_stats.csv",
            category="lineage"),
        "{outdir}/enrichment/shuffle_stats/histogram_p_values_max.png"
    params:
        outdir = lambda wildcards, output: dirname(output[0]),
        note = join("src", "clones", "wrap_clones_enrich_shuffle.ipynb"),
        samples=",".join(config["samples"].index),
    threads: 16
    shell: "papermill -p enrich_f {input.enrich_f} -p cells_meta_f {input.cells_meta_f} -p OUTDIR {params.outdir} -p samples {params.samples:q} {params.note} {output[0]}"


rule finalize:
    input:
        "{outdir}/enrichment/enrichmentNorm.csv",
        "{outdir}/enrichment/shuffle_stats/shuffle_stats.csv",
    output:
        "{outdir}/enrichment/_enrichment_complete",
    shell: "touch {output}"