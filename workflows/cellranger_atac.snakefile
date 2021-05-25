####
## WILL CHANGE TO THE DIRECTORY listed as mtscATAC_OUTDIR or OUTDIR in configfile.
## Directory will be made if not already there. Will load the samples csv file before
## changing directory.
####
import os
from os.path import abspath, dirname
from shutil import rmtree
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyBigWig
mpl.use('Agg')
import numpy as np
samples_df = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["out_name"], drop=False)
samples_df.columns = samples_df.columns.str.strip()
print(samples_df)

if "mtscATAC_OUTDIR" in config:
    outdir = config["mtscATAC_OUTDIR"]
elif "OUTDIR" in config:
    outdir = config["OUTDIR"]
else:
    raise ValueError("OUTDIR not in config")
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.chdir(outdir)
##
#report: "report/workflow.rst"
rule all:
    input:
        expand("{sample}/outs/possorted_bam.bam", sample=samples_df.index),
        expand("{sample}/outs/web_summary.html", sample=samples_df.index),
        #"aggregate/outs/web_summary.html",
        "reanalysis/outs/web_summary.html",
        expand("coverage/{sample}_coverage.bw", sample=samples_df.index),
         #expand("coverage/{sample}_coverage.tsv", sample=samples_df.index),
        expand("plots/{sample}_coverage_chr.png", sample=samples_df.index),


def results_dir():
    return config["results"]


rule cellranger_count:
    params:
        fastq_path = config['fastq_dir'],  #lambda wildcards: samples_df.loc[wildcards.out_name, 'fastq_dir'],
        sample = lambda wildcards: samples_df.loc[wildcards.sample, 'sample'],
        out_name = lambda wildcards: wildcards.sample,
        genome_ref = config["genome_dir"]
    output:
        bam_f = "{sample}/outs/possorted_bam.bam",
        html_out = "{sample}/outs/web_summary.html",
        frag_out = "{sample}/outs/fragments.tsv.gz",
        sing_out = "{sample}/outs/singlecell.csv"
    threads: 25
    #log: "logs/cellratac_{sample}.log "
    run:
        f = os.path.dirname(output.bam_f).split("/")[0]
        print('f', f)
        print('cwd', os.getcwd())
        rmtree(f)
        print(os.getcwd())
        print(os.listdir('.'))
        shell("cellranger-atac count --id={params.out_name} --localmem=115 --reference={params.genome_ref} --fastqs={params.fastq_path} --sample={params.sample} --localcores={threads}")


rule generate_coverage:
    input: "{sample}/outs/possorted_bam.bam"
    output:
        coverage="coverage/{sample}_coverage.tsv",
    shell: "samtools index {input};samtools depth {input} > {output} "


rule generate_bigwig_coverage:
    input: "{sample}/outs/possorted_bam.bam"
    output:
        coverage="coverage/{sample}_coverage.bw"
    shell: "samtools index {input};bamCoverage -b {input} -o {output}"

# rule chrom_coverage:
#     input:
#         cov_in = "coverage/{sample}_coverage.tsv"
#     output:
#         chrom_cov = report("plots/{sample}_coverage_chr.png", caption="report/fig1.rst", category="Coverage")
#     run:
#         print(input.cov_in)
#         cov = pd.read_csv(input.cov_in, sep='\t', header=None)
#         cov.columns = ["Chr", "Pos", "Cov"]
#         # If GRCh38, drop any chr with "_" in it,
#         # which are scaffolds.
#         cov = cov[~(cov["Chr"].str.contains("_"))]
#         cov_chr = cov.groupby("Chr").sum()["Cov"].reset_index()
#         cov_chr.to_csv(output.chrom_cov+".csv")
#         cov_chr["log10 coverage"] = np.log10(cov_chr["Cov"]+1)
#         sns.barplot(x="Chr", y="log10 coverage", data=cov_chr)
#         plt.xlabel("Chromosome")
#         plt.xticks(rotation=90)
#         plt.savefig(output.chrom_cov)
#         plt.close()


rule chrom_coverage:
    "Average coverage for chromosome"
    input: rules.generate_bigwig_coverage.output
    output:
        chrom_cov = report("plots/{sample}_coverage_chr.png", caption="report/fig1.rst", category="Coverage")
    run:
        all_chroms={}
        bw = pyBigWig.open(input[0])
        for chrom in bw.chroms():
             all_chroms[chrom] = bw.stats(chrom)[0]
        all_chroms = pd.Series(all_chroms , name='Chromosome average reads per bp')
        all_chroms=all_chroms.loc[~(all_chroms.index.str.contains('chrUn'))]
        sns.barplot(x=all_chroms.index,y=all_chroms.values)
        plt.xlabel("Chromosome")
        plt.xticks(rotation=90)
        plt.savefig(output[0])
        plt.close()
        sns.barplot(x=all_chroms.index,y=np.log10(all_chroms.values+1))
        plt.xlabel("Chromosome")
        plt.xticks(rotation=90)
        plt.savefig(output[0].replace('.png','')+'.log10.png')
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


rule create_csv:
    input:
        expand("{sample}/outs/fragments.tsv.gz", sample=samples_df.index),
    output:
        aggr_csv = "aggr.csv"
    run:
        df = pd.DataFrame(columns=["library_id", "fragments", "cells"])
        print('in',input)
        for ind, frag in enumerate(input):
            print('frag', frag)
            cell = frag.replace("fragments.tsv.gz", "singlecell.csv")
            samp = os.path.dirname(frag).split('/')[0]
            print(pd.DataFrame([samp, frag, cell], index=["library_id", "fragments", "cells"]).transpose())
            df = pd.concat((df,
                           pd.DataFrame([samp, os.path.abspath(frag), os.path.abspath(cell)], index=["library_id", "fragments", "cells"]).transpose()),
                           axis=0)
            print('df', df)
        df.to_csv(output.aggr_csv, index=False)

#        library_id,fragments,cells
#LV123,/opt/runs/LV123/outs/fragments.tsv.gz,/opt/runs/LV123/outs/singlecell.csv
#LB456,/opt/runs/LB456/outs/fragments.tsv.gz,/opt/runs/LB456/outs/singlecell.csv
    #params: lambda wildcards: wildcards.sa

rule aggr_cellranger:
    input:
        aggr_csv = "aggr.csv" #expand("{sample}/outs/web_summary.html", out_)
    output: "aggregate/outs/web_summary.html"
    params:
        ref = config["genome_dir"],
    shell: "cellranger-atac aggr --id=aggregate --csv={input.aggr_csv} --localmem=50 --normalize=depth --reference={params.ref}"

rule reanalyze:
    input: "aggregate/outs/web_summary.html"
    output: "reanalysis/outs/web_summary.html"
    params:
        aggr_outs = "aggregate/outs",
        ana_outs = "reanalysis/outs",
        ref = config["genome_dir"],
        full_in = lambda wildcards, input: dirname(abspath(input[0])),
        full_out = lambda wildcards, output: dirname(abspath(output[0]))
    threads: 24
    shell:
         "cellranger-atac reanalyze --id=reanalysis \
          --peaks={params.full_in}/peaks.bed \
          --reference={params.ref} --localcores={threads} --localmem=50 \
          --fragments={params.full_in}/fragments.tsv.gz"
          #--params={params.full_out}/aggregation_csv.csv \