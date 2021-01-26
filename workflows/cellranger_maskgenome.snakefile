import os
os.chdir(config["OUTDIR"])

rule all:
    input:
        config["name"],
        ".snake_completed",


rule mask_mt:
    params:
        fasta = config["genome_in"],
        blacklist = config["blacklist_f"]
    output: config["genome_blacklist_f"]
    shell: "bedtools maskfasta -fi {params.fasta} -bed {params.blacklist} -fo {output}"


rule cellranger_mkref:
    input: rules.mask_mt.output
    params:
        name = config["name"],
        config = config["config_file"]
    output:
        directory(config["name"]),
        ".snake_completed"
    run:
         shell("cellranger-atac mkref {params.name} --config {params.config}"),
         shell("touch {output}")

