wildcard_constraints:
    mapq="\d+",
    cellr='True|False'

from src.config import r_kernel
import os
import pandas as pd
from snakemake.utils import validate
from Bio import SeqIO
import pickle
import subprocess as subp
from datetime import datetime
import copy

# Flatten dictionary parameters
if "mttrace" in config:
    for c in config['mttrace']:
        config[c] = copy.deepcopy(config['mttrace'][c])
# Flatten dictionary parameters
if "mgatk" in config:
    for c in config['mgatk']:
        config[c] = copy.deepcopy(config['mgatk'][c])

print('config')
print(config)
cellr_bc = config["use_cellr_barcode"]
num_reads_filter = config["num_reads_filter"]
maxBP = config["maxBP"]
ref_fa = config["ref_fa"]
#print(pd.read_table(config["samples"], dtype=str,sep=','))
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)


#print('index',samples.index)
res = config["results"]
#mq = config["mapq"]

ft = config["filters"]
#workdir: config["work_dir"]
rule all:
    input:
        expand("{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
               results=res,sample=samples["sample_name"].values, num_read=config["num_reads_filter"]),
        # expand("{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/{sample}.lowC{low_cov_thresh}.variant.rds",
        #        results=res,sample=samples["sample_name"].values,
        #        cellr_bc=cellr_bc, num_read=num_reads_filter, low_cov_thresh=config["low_cov_thresh"]),
        expand("{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png",
               results=res,sample=samples["sample_name"].values, num_read=config["num_reads_filter"], cellr_bc=config["use_cellr_barcode"])



def get_sample_bam(wildcards):
    bam = samples.loc[wildcards.sample,"bam_f"]
    #print(bam)
    return bam # + ".bam"


def get_sample_barcodes(wildcards):
    return samples.loc[wildcards.sample, "barcode_f"]
    #return os.path.join(os.path.dirname(samples.loc[wildcards.sample,"bam_f"]), "filtered_feature_bc_matrix/barcodes.tsv")


def results_dir():
    return config["results"]


rule link_bam:
    input: get_sample_bam
    output: "{results}/{sample}/00_bam/{sample}.bam"
    shell: "ln -s {input} {output}"

rule index_bam:
    """Index the bam file"""
    input: rules.link_bam.output
    output: "{results}/{sample}/00_bam/{sample}.bam.bai"
    shell: "samtools index {input}"


# rule filter_bam:
#     """Filter reads with low MAPQ and other potential filters """
#     input:
#         get_sample_bam,
#         #in_bam = "{results}/{sample}/00_bam/{sample}.bam",
#         #in_bai = "{results}/{sample}/00_bam/{sample}.bam.bai",
#     output:
#         mq_bam = "{results}/{sample}/01_bam_filter/{sample}.bam",
#         mq_bai = "{results}/{sample}/01_bam_filter/{sample}.bam.bai"
#
#     run:
#             if "rm_duplicates" in config and config["rm_duplicates"]:
#                 shell("samtools view {input.in_bam} -F 1284 -q {wildcards.mapq} -h -b > {output.mq_bam}")
#             else:
#                 shell("samtools view {input.in_bam} -q {wildcards.mapq} -h -b > {output.mq_bam}")
#             shell("samtools index {output.mq_bam}")


rule MT_map:
    """Extract the MT genome"""
    # input:
    #     bam = "{results}/{sample}/01_bam_filter/{sample}.bam",
    #     bai = "{results}/{sample}/01_bam_filter/{sample}.bam.bai"
    input:
        #bam = "{results}/{sample}/00_bam/{sample}.bam",
        bai = "{results}/{sample}/00_bam/{sample}.bam.bai",
    output:
        mt_bam="{results}/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/{sample}/MT/{sample}.MT.bam.bai"
    params:
        bam = lambda wildcards, input: input.bai.replace('.bai', ''),
        mt_chr=config["mito_character"],

    run:
        #shell("samtools {input.bam}")
        shell("samtools view -b {params.bam} {params.mt_chr} > {output.mt_bam}")
        shell("samtools index {output.mt_bam}")

rule get_bigwig:
    """Extract the MT genome"""
    input: "{results}/{sample}/MT/{sample}.MT.bam"
    output:
        coverage="{results}/{sample}/MT/{sample}.MT.bw"
    shell: "bamCoverage -b {input} -o {output}"


rule barcode_data:
    """Loop through the bam file and extract the barcode information."""
    input:
        mt_bam="{results}/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/{sample}/MT/{sample}.MT.bam.bai"
    output: "{results}/{sample}/MT/{sample}_barcode_data.p",
    params: mt_chr=config["mito_character"]
    shell:
         "python src/bam_barcodes_function.py {input.mt_bam} {output} {params.mt_chr}"


rule barcode_filter:
    input:
        barcode_p = "{results}/{sample}/MT/{sample}_barcode_data.p",
        cellr_f = get_sample_barcodes
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    params: cellr_bc = "{cellr_bc}"
    shell: "python src/filter_barcodes.py {input} {params} {output}"


rule sortCB:
    input:
        mt_bam="{results}/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/{sample}/MT/{sample}.MT.bam.bai"
    output: "{results}/{sample}/MT/{sample}.MT.CB.bam"
    shell: "samtools sort -t CB {input.mt_bam} > {output}"


rule scBam:
    """Extract each single-cell and put into respective bam file"""
    input: "{results}/{sample}/MT/{sample}.MT.CB.bam"
        #mt_bam="{results}/{sample}/MT/{sample}.MT.bam",
        #mt_bai="{results}/{sample}/MT/{sample}.MT.bam.bai"
    output: directory("/data/isshamie/mito_lineage/{results}/{sample}/MT/{sample}_scBam")
    threads: 18
    #directory("/data/isshamie/miito_lineage/mttrace/{sample}/MT/{sample}_scBam")
    shell:
        "python src/split_by_CB.py {input} {output}"


rule scPileup:
    """Run the first part of the MT-genotype function by getting the read pileups for each bam file for each nucleotide and overall coverage"""
    input:
        scBam = "/data/isshamie/mito_lineage/{results}/{sample}/MT/{sample}_scBam",
        #scBam = "{results}/{sample}/MT/{sample}_scBam",
        barcodes = "{results}/{sample}/MT/{sample}_barcode_data.p",
    output:
          directory("{results}/{sample}/MT/{sample}_scPileup_{num_read}")
    params:
        base_quality = config['base_quality']
    shell:
         "python src/scPileup_counts.py {input.scBam} {output} {input.barcodes} {wildcards.num_read} {params.base_quality}"

def concat_files(directory, samplename, nt):
    cmd = f"find {directory} -type f -name *.{nt}.txt -exec cat {{}} \; > {samplename}_all.{nt}.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)
    cmd = f"find {directory} -type f -name *.{nt}.minus.txt -exec cat {{}} \; > {samplename}_all.{nt}.minus.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)
    return

rule scConcat:
    input:
        scPileup_dir = "{results}/{sample}/MT/{sample}_scPileup_{num_read}",
        scBam = "/data/isshamie/mito_lineage/{results}/{sample}/MT/{sample}_scBam"
    output:
        "{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.minus.txt"
    threads: 28
    params:
        samplename = lambda wildcards, output: output[0].split("_all.coverage.minus.txt")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    run:
        for n in ["A", "C", "G", "T", "coverage"]:
            concat_files(input.scPileup_dir, params.samplename, n)


rule scPileup_concat_strands:
    """ Run the second part of the MT-genotype pipeline, which just concatenates all the pileup data for each nucleotide and overall."""
    input:  "{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.minus.txt"
    output:
        all = "{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz"
    params:
        concat_dir = lambda wildcards, input: os.path.dirname(input[0]),
        samplename = lambda wildcards, output: output.all.split("_all.coverage.strands.txt.gz")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    shell:
         "python src/scPileup_concat.py {params.concat_dir} {params.samplename}"

# scPileup_concat with different strands
# rule scPileup_concat_strands:
#     """ Run the second part of the MT-genotype pipeline, which just concatenates all the pileup data for each nucleotide and overall."""
#     input:
#         scPileup_dir = "{results}/{sample}/MT/{sample}_scPileup_{num_read}",
#         # Added to save directly to /data/isshamie
#         scBam = "/data/isshamie/mito_lineage/{results}/{sample}/MT/{sample}_scBam"
#     output:
#         all = "{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz"
#     threads: 34
#     params:
#         samplename = lambda wildcards, output: output.all.split("_all.coverage.strands.txt.gz")[0]
#           #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
#     shell:
#          "python src/scPileup_concat.py {input.scPileup_dir} {params.samplename}"

rule plot_CB_coverage:
    """Plot the MT coverage of the single cells"""
    input:
         all = "{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
         barcodes = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_CB_coverage_hist_minReads{num_read}.png"
    shell:
        "python src/plot_CB_coverage.py {input} {output}"


rule scPileup_MT_matrix:
    """Create the position-by-cell coverage matrix"""
    input:
        all = "{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
        barcode_p = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output:
        sc_coverage_f = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    threads: 16
    shell:
        "python src/plot_heatmap_coverage.py sc_mt {input.barcode_p} {input.all} {output.sc_coverage_f} {maxBP}"


def run_filter_cell(in_f, prefix, barcode_p, n=""):
    df = pd.read_csv(in_f)
    df["CB"] = df["CB"].str.replace(".bam","")
    barcodes = pickle.load(open(barcode_p, "rb"))
    if isinstance(barcodes, dict):
        df = df[df["CB"].isin(list(barcodes.keys()))]
    else:
        df = df[df["CB"].isin(barcodes)]
    if n == "coverage":
        df = df.iloc[:,:3]
    df.to_csv(prefix+f".{n}.strands.txt.gz", header=None, index=None, compression='gzip')
    return


rule filter_cell_bc:
    """Extracts only the relevant cell barcodes and removes the .bam from the barcode names."""
    input:
        all = "{results}/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
        barcode_p = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output:
        "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz"
    run:
        for n in ["A", "C", "G", "T", "coverage"]:
            curr_f = input.all
            curr_f = curr_f.replace(".coverage.", "." + n + ".")
            df = pd.read_csv(curr_f)
            df["CB"] = df["CB"].str.replace(".bam","")
            barcodes = pickle.load(open(input.barcode_p, "rb"))
            if isinstance(barcodes, dict):
                df = df[df["CB"].isin(list(barcodes.keys()))]
            else:
                df = df[df["CB"].isin(barcodes)]
            curr_out_f = output[0]
            curr_out_f = curr_out_f.replace(".coverage.", "." + n + ".")
            if n == "coverage":
                df = df.iloc[:,:3]
            df = df.sort_values(["CB", "Position"])
            df.to_csv(curr_out_f, header=None, index=None, compression='gzip')


def get_ref(wildcards):
    w = wildcards
    mito = config['mito_character']
    return f"{w.results}/{w.sample}/mapq_{w.mapq}/cellr_{w.cellr_bc}/{w.sample}_{w.num_read}/{mito}_refAllele.txt"


rule plot_scPileup_MT_matrix:
    """Plot the posiitonal coverages."""
    input:
        sc_coverage_f = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    output:
        save_f_heat = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position.png",
        save_f_coverage = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png"
    shell:
        #"python src/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_coverage}"
        "python src/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_heat} {output.save_f_coverage}"


def get_filt(w):
    return w.min_cells, w.min_reads, w.topN, w.het_thresh, w.min_het_cells, w.het_count_thresh, w.bq_thresh


rule create_filters:
    input:
        concat_dir = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz"
    #output:  "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/af_by_cell.tsv"
    output:
        af_f = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/af_by_cell.tsv",
        cov = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt"
    params:
        concat_d = lambda wildcards, input: os.path.dirname(input.concat_dir),
        ref_fa = config['mt_ref_fa'],
        name = lambda wildcards: wildcards.sample,
        filt_params = get_filt
    resources:
        mem_mb=90000
    log: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/.outcfg"
    shell: "python src/calculate_AF_by_cell.py {params.concat_d} {output.af_f} {params.ref_fa} {params.name} {params.filt_params} --log {log}"


rule get_refAllele:
    #input: config["mt_ref_fa"],
    params: config["chrM_refAllele"]
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/chrM_refAllele.txt"
    shell: 'cp {input} {output}'


rule to_seurat:
    """ Extract the variants based on the specified parameters"""
    input:
        all = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
        #all = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz",
        #depth = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.depthTable.txt",
        refAllele = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/chrM_refAllele.txt"
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.rds"
    params:
        data_dir=lambda wildcards, input: os.path.dirname(input.all),
        sample = lambda wildcards: wildcards.sample,
        rkernel = r_kernel
    log: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/.to_seurat.log"
    shell:
        "./R_scripts/toRDS {params.data_dir} filter_mgatk/{params.sample} FALSE" # &> {log}"


rule callVariants_mgatk:
    input: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.rds"
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/{sample}.lowC{low_cov_thresh}.variant.rds"
    params:
        low_cov_thresh=lambda wildcards: wildcards.low_cov_thresh,
        rkernel = r_kernel
    #log: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/mgatk/.callVariants_mgatk.lowC{low_cov_thresh}.log"
    shell: "./R_scripts/variant_calling {input} {params.low_cov_thresh}"# &> {log}"



# B. Run mgatk before the filters
rule raw_get_depth:
    input:  "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz" #"{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.depthTable.txt"
    run:
        df = pd.read_csv(input[0], header=None)
        depth = df.groupby(1).sum()[2]/16569
        depth.to_csv(output[0], sep='\t', header=False)


rule raw_get_refAllele:
    input: config["mt_ref_fa"],
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/chrM_refAllele.txt"
    run:
        records = list(SeqIO.parse(input[0], "fasta"))
        mt_seq = records[0].seq
        for ind, val in enumerate(mt_seq):
            if ind == 0:
                ref_str = f"{ind+1}\t{val}"
            else:
                ref_str = f"{ref_str}\n{ind+1}\t{val}"
        with open(output[0], "w") as f:
            f.write(ref_str)

rule raw_to_seurat:
    """ Extract the variants based on the specified parameters"""
    input:
        #all = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.txt",
        all = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz",
        depth = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.depthTable.txt",
        refAllele = "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/chrM_refAllele.txt"
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/{sample}.rds"
    params:
        data_dir=lambda wildcards, input: os.path.dirname(input.all),
        sample = lambda wildcards: wildcards.sample,
        rkernel = r_kernel
    log: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/.to_seurat.log"
    shell:
        "./R_scripts/toRDS {params.data_dir} mgatk/{params.sample} TRUE" # &> {log}"


rule raw_callVariants_mgatk:
    input:  "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/{sample}.rds"
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/{sample}.lowC{low_cov_thresh}.variant.rds"
    params:
        rkernel = r_kernel,
        low_cov_thresh=lambda wildcards: wildcards.low_cov_thresh
    log: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/mgatk/.callVariants_mgatk.lowC{low_cov_thresh}.log"
    shell: "./R_scripts/variant_calling {input} {params.low_cov_thresh}" # &> {log}"



rule plot_cells:
    input: "{results}/cellr_{cellr_bc}_nr_{num_reads}_mc_{min_cells}_mr_{min_reads}/aggregate_af.csv"
    output: "{results}/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/cluster_{min_cells}_{min_reads}_{cell_mt_coverage}.png"
    shell: "python plot_lineage.py {input} {output}"



