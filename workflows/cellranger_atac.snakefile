####
## WILL CHANGE TO THE DIRECTORY listed as OUTDIR in configfile.
## Directory must already be made. Will load the samples csv file before
## changing directory.
####
import os
from shutil import rmtree
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
samples_df = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["out_name"], drop=False)
samples_df.columns = samples_df.columns.str.strip()
print(samples_df)
os.chdir(config["OUTDIR"])
##


#report: "report/workflow.rst"


rule all:
    input:
        expand("{sample}/outs/possorted_bam.bam", sample=samples_df.index),
        expand("{sample}/outs/web_summary.html", sample=samples_df.index),
        expand("coverage/{sample}_coverage.tsv", sample=samples_df.index),
        expand("plots/{sample}_coverage_chr.png", sample=samples_df.index)

def results_dir():
    return config["results"]


rule cellranger_count:
    params:
        fastq_path = lambda wildcards: samples_df.loc[wildcards.out_name, 'fastq_dir'],
        sample = lambda wildcards: samples_df.loc[wildcards.out_name, 'sample'],
        out_name = lambda wildcards: wildcards.out_name,
        genome_ref = config["genome_dir"]

    output:
        bam_f = "{out_name}/outs/possorted_bam.bam",
        html_out = "{out_name}/outs/web_summary.html"
    threads: 18
    run:
        f = os.path.dirname(output.bam_f).split("/")[0]
        print('f', f)
        print('cwd', os.getcwd())
        rmtree(f)
        print(os.listdir('.'))
        shell("cellranger-atac count --id={params.out_name} --localmem=115 --reference={params.genome_ref} --fastqs={params.fastq_path} --sample={params.sample} --localcores={threads}")


rule generate_coverage:
    input: "{out_name}/outs/possorted_bam.bam"
    output:
        coverage="coverage/{out_name}_coverage.tsv",
    shell: "samtools depth {input} > {output} "


rule chrom_coverage:
    input:
        cov_in = "coverage/{out_name}_coverage.tsv"
    output:
        chrom_cov = report("plots/{out_name}_coverage_chr.png", caption="report/fig1.rst", category="Coverage")
    run:
        print(input.cov_in)
        cov = pd.read_csv(input.cov_in, sep='\t', header=None)
        cov.columns = ["Chr", "Pos", "Cov"]
        # If GRCh38, drop any chr with "_" in it,
        # which are scaffolds.
        cov = cov[~(cov["Chr"].str.contains("_"))]
        cov_chr = cov.groupby("Chr").sum()["Cov"].reset_index()
        cov_chr.to_csv(output.chrom_cov+".csv")
        cov_chr["log10 coverage"] = np.log10(cov_chr["Cov"]+1)
        sns.barplot(x="Chr", y="log10 coverage", data=cov_chr)
        plt.xlabel("Chromosome")
        plt.xticks(rotation=90)
        plt.savefig(output.chrom_cov)
        plt.close()


rule cellranger_qc:
    input: expand("{sample}/outs/web_summary.html", sample=samples_df.index)
    output: "plots/cellr_qc.png"
    shell: "python qc_cellr.py {input}"
