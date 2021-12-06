################################
## Comparing algorithmically called clones
################################
rule aggregate_clones_qc:
    input: expand("{{outdir}}/{cluster}/{cluster_params}/cells_meta.tsv",
                   cluster=config["cluster"]["params"],
                   cluster_params=config["cluster"]["params"] )

rule compare_clones:
    input: expand("{{outdir}}/{cluster}/{cluster_params}/cells_meta.tsv",
               cluster=config["cluster"]["params"],
               cluster_params=config["cluster"]["params"])
    run:
        from src.clones.compare_clones import compare_clones
        compare_clones()

################################
## Comparing algorithmically called variants
################################
rule compare_variants:
    input: expand("{{outdir}}/{cluster}/{cluster_params}/cells_meta.tsv",
                   cluster=config["cluster"]["params"],
                   cluster_params=config["cluster"]["params"] )

