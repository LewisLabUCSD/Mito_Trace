configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"

rule all:
    input:

rule orig_index:
    input: config["raw_folder"] / config["bam_f"]
         #{raw_folder}/{bam_f}"
    output: "data/processed/{name}.bam.bai"
    shell: "samtools index {input}"

rule MT_map:
    input:
         "{raw_folder}/{bam_f}"
    output:
          bai="data/processed/{name}.bam.bai"
          mt_bam="data/processed/{name}.MT.bam"
          mt_bai="data/processed/{name}.MT.bam.bai"
    shell:
         "samtools view -b {input} MT > {output.mt_bam}"
         "samtools index {output.mt_bam}"

rule barcode_data:
    input: "data/processed/{name}.MT.bam"
    output: "data/processed/{name}_barcode_data.p"
    shell:
         "python src/bam_barcodes_function.py {input} {output}"

#rule: