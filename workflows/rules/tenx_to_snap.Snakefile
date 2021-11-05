import os
#os.chdir(config["OUTDIR"])

genome_name = config["genome"]
genome_size = config["genome_size"]


rule all:
    input: expand("data/{name}_cell_by_bin.txt", name=config["samples"])

def get_bam(wildcards):
    return os.path.join(config["indir"], wildcards.name, "outs/possorted_bam.bam")


rule sortname_bam:
    input: get_bam
    output: "data/{name}.nsrt.bam"
    shell: "samtools sort -n -@ 10 -m 1G {input} -o {output}"

rule snap_pre:
    input:  "data/{name}.nsrt.bam"
    output:
        #snap="{name}.snap", This is always overwritten, so not great for file checking
        complete ="{name}_snap_pre.txt"
    params:
        name = lambda wildcards: wildcards.name,
        genome_size = genome_size,
        genome_name = genome_name
    resources:
        mem_mb=100000
    log: "logs/snap_pre_{name}.log"
    run:
        shell("snaptools snap-pre  \
	       --input-file={input}    \
	       --output-snap=data/{params.name}.snap    \
	       --genome-name={params.genome_name}  \
	       --genome-size={params.genome_size}  \
	       --min-mapq=30  \
	       --min-flen=0  \
	       --max-flen=1000  \
	       --keep-chrm=TRUE  \
	       --keep-single=TRUE  \
	       --keep-secondary=False  \
	       --overwrite=True  \
	       --max-num=1000000  \
	       --min-cov=100  \
	       --verbose=True 2> {log}"),
        shell("touch {output.complete}")


rule snap_add_cell_by_bin:
    input:
        complete ="{name}_snap_pre.txt"
    output:
        #snap="{name}.snap", This is always overwritten, so not great for file checking
        complete =  "data/{name}_cell_by_bin.txt"
    params:
        name = lambda wildcards: wildcards.name
    log: "logs/snap_add_cell_by_bin_{name}.log"
    resources:
        mem_mb=100000
    run:
        shell("snaptools snap-add-bmat  \
	        --snap-file=data/{params.name}.snap \
	        --bin-size-list 5000 10000  \
	        --verbose=True 2> {log}"),
        shell("touch {output.complete}")
