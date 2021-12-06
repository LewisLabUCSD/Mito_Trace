
import pandas as pd


samples = pd.read_table(config["samples"], sep=',')
samples.columns = samples.columns.str.strip()
samples["sra_id"] = samples["sra_id"].astype(str).str.strip()
samples["name"] = samples["name"].astype(str).str.strip()
sample_names = samples["name"]
samples = samples.set_index("name")
#print(samples)
n_cores = config["n_cpu"]
max_mem = config["max_mem_kb"]
sra_id = samples["sra_id"]
print(samples)

report: "report/workflow.rst"

rule all: 
	input: 
		expand("out/{sample_name}_1.fastq.gz", sample_name=sample_names),
		expand("out/{sample_name}_2.fastq.gz", sample_name=sample_names),
		expand("out/{sample_name}_3.fastq.gz", sample_name=sample_names)

rule mk_temp:
	output: directory("temp")
	shell: "mkdir {output}"

def get_sra(wildcards):
	print('wildcards')
	print(wildcards)
	return samples.loc[wildcards.sample_name, "sra_id"]


def get_sra_in(wildcards):
	sra_id = samples.loc[wildcards.sample_name, "sra_id"]
	print(sra_id)
	return f"{sra_id}_1.fastq.gz", f"{sra_id}_2.fastq.gz", f"{sra_id}_3.fastq.gz"


rule download_sra_to_fastq:
	input: 
		rules.mk_temp.output
	output:
		"out/{sample_name}_1.fastq.gz",
		"out/{sample_name}_2.fastq.gz",
		"out/{sample_name}_3.fastq.gz" 
	threads: n_cores
	params: 
		sra=get_sra,
		sra_in = get_sra_in,
		max_mem=max_mem
	run: 
		shell("prefetch {params.sra} --max-size {params.max_mem}")
		shell("parallel-fastq-dump --sra-id {params.sra} --threads {threads} --gzip --split-files --tmpdir {input}")	
		shell("mv {params.sra_in[0]} {output[0]}"),
		shell("mv {params.sra_in[1]} {output[1]}"),
		shell("mv {params.sra_in[2]} {output[2]}")

