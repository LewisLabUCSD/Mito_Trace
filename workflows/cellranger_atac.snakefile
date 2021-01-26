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
if not os.path.exists(config["OUTDIR"]):
    os.mkdir(config["OUTDIR"])
os.chdir(config["OUTDIR"])
##


#report: "report/workflow.rst"


rule all:
    input:
        expand("{sample}/outs/possorted_bam.bam", sample=samples_df.index),
        expand("{sample}/outs/web_summary.html", sample=samples_df.index),
        expand("coverage/{sample}_coverage.tsv", sample=samples_df.index),
        expand("plots/{sample}_coverage_chr.png", sample=samples_df.index),
        expand("coverage/{sample}_coverage.bw", sample=samples_df.index)

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
    threads: 24
    #log: "logs/cellratac_{out_name}.log "
    run:
        f = os.path.dirname(output.bam_f).split("/")[0]
        print('f', f)
        print('cwd', os.getcwd())
        rmtree(f)
        print(os.getcwd())
        print(os.listdir('.'))
        shell("cellranger-atac count --id={params.out_name} --localmem=115 --reference={params.genome_ref} --fastqs={params.fastq_path} --sample={params.sample} --localcores={threads}")


rule generate_coverage:
    input: "{out_name}/outs/possorted_bam.bam"
    output:
        coverage="coverage/{out_name}_coverage.tsv",
    shell: "samtools index {input};samtools depth {input} > {output} "


rule generate_bigwig_coverage:
    input: "{sample}/outs/possorted_bam.bam"
    output:
        coverage="coverage/{sample}_coverage.bw"
    shell: "samtools index {input};bamCoverage -b {input} -o {output}"

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

rule flagstat:
    input: "{sample}/outs/possorted_bam.bam"
    output: "{sample}/outs/possorted_bam.flagstat.txt"
    shell: "samtools flagstat {input} > {output}"

rule mapq:
    input: "{sample}/outs/possorted_bam.bam"
    output: "{sample}/outs/possorted_bam.mapq.txt"
    shell: """samtools view {input} | awk -F "\t" '{print $5}' > {output} """

rule cellranger_qc:
    """TODO
    """
    input: expand("{sample}/outs/web_summary.html", sample=samples_df.index)
    output: "plots/cellr_qc.png"
    shell: "python qc_cellr.py {input}"


rule subsample_fastq:
    shell: """paste {input} | awk '{ printf("%s",$0); n++;
    if(n%4==0) { printf("\n");} else { printf("\t");} }' |
    awk -v k=10000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |
    awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "forward_sub.fastq";print
    $2"\n"$4"\n"$6"\n"$8 > "reverse_sub.fastq"}'
    """
