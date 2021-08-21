import pandas as pd
import pickle
from os.path import dirname
from snakemake.utils import min_version
min_version("6.0")

res = config['results']
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)


mt = config['mtpreproc']
ft = config['filters']
rule all:
    input:
        expand("{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
            results=res,
            sample=samples["sample_name"].values,
            num_read=config["num_reads_filter"]),
        expand("{results}/data/{sample}/MT/cellr_{cellr_bc}/numread_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt",
               results=res, sample=samples["sample_name"].values,
               num_read=mt['num_reads_filter'], cellr_bc=mt['cellr_bc'],
               min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']['n_clone_list'])


def get_sample_bam(wildcards):
    bam = samples.loc[wildcards.sample,"bam_f"]
    print(bam)
    return bam


def get_sample_barcodes(wildcards):
    return samples.loc[wildcards.sample, "barcode_f"]


rule link_bam:
    input: get_sample_bam
    output: "{results}/data/{sample}/00_bam/{sample}.bam"
    shell: 'ln -sr {input} {output}'

rule index_bam:
    """Index the bam file"""
    input: "{results}/data/{sample}/00_bam/{sample}.bam"#rules.link_bam.output
    output: "{results}/data/{sample}/00_bam/{sample}.bam.bai"
    shell: "samtools index {input}"


rule MT_map:
    """Extract the MT genome"""
    input:
        #bam = "{results}/data/{sample}/00_bam/{sample}.bam",
        bai = "{results}/data/{sample}/00_bam/{sample}.bam.bai",
    output:
        mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    params:
        bam = lambda wildcards, input: input.bai.replace('.bai', ''),
        mt_chr=config["mito_character"],

    run:
        shell("samtools view -b {params.bam} {params.mt_chr} > {output.mt_bam}")
        shell("samtools index {output.mt_bam}")

rule get_bigwig:
    """Extract the MT genome"""
    input: "{results}/data/{sample}/MT/{sample}.MT.bam"
    output:
        coverage="{results}/data/{sample}/MT/{sample}.MT.bw"
    shell: "bamCoverage -b {input} -o {output}"


rule barcode_data:
    """Loop through the bam file and extract the barcode information."""
    input:
        mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: "{results}/data/{sample}/MT/{sample}_barcode_data.p",
    params: mt_chr=config["mito_character"]
    shell:
         "python src/bam_barcodes_function.py {input.mt_bam} {output} {params.mt_chr}"


rule barcode_filter:
    input:
        barcode_p = "{results}/data/{sample}/MT/{sample}_barcode_data.p",
        cellr_f = get_sample_barcodes
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    params: cellr_bc =  lambda wildcards: wildcards.cellr_bc
    shell: "python src/filter_barcodes.py {input} {params} {output}"


rule sortCB:
    input:
        mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: "{results}/data/{sample}/MT/{sample}.MT.CB.bam"
    shell: "samtools sort -t CB {input.mt_bam} > {output}"


rule scBam:
    """Extract each single-cell and put into respective bam file"""
    input: "{results}/data/{sample}/MT/{sample}.MT.CB.bam"
        #mt_bam="{results}/data/{sample}/MT/{sample}.MT.bam",
        #mt_bai="{results}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: directory("{results}/data/{sample}/MT/{sample}_scBam")
    threads: 18
    #directory("/data/isshamie/miito_lineage/mttrace/{sample}/MT/{sample}_scBam")
    shell:
        "python src/split_by_CB.py {input} {output}"


rule scPileup:
    """Run the first part of the MT-genotype function by getting the read pileups for each bam file for each nucleotide and overall coverage"""
    input:
        scBam = "{results}/data/{sample}/MT/{sample}_scBam",
        #scBam = "{results}/data/{sample}/MT/{sample}_scBam",
        barcodes = "{results}/data/{sample}/MT/{sample}_barcode_data.p",
    output:
          directory("{results}/data/{sample}/MT/{sample}_scPileup_{num_read}")
    params:
        base_quality = config['base_quality']
    threads: 18
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
        scPileup_dir = "{results}/data/{sample}/MT/{sample}_scPileup_{num_read}",
        scBam = "{results}/data/{sample}/MT/{sample}_scBam"
    output:
        "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.minus.txt"
    threads: 28
    params:
        samplename = lambda wildcards, output: output[0].split("_all.coverage.minus.txt")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    run:
        for n in ["A", "C", "G", "T", "coverage"]:
            concat_files(input.scPileup_dir, params.samplename, n)


rule scPileup_concat_strands:
    """ Run the second part of the MT-genotype pipeline, which just concatenates all the pileup data for each nucleotide and overall."""
    input:  "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.minus.txt"
    output:
        all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz"
    params:
        concat_dir = lambda wildcards, input: dirname(input[0]),
        samplename = lambda wildcards, output: output.all.split("_all.coverage.strands.txt.gz")[0]
          #"data/processed/{sample}/scPileup_concat/{sample}_{num_read}"
    shell:
         "python src/scPileup_concat.py {params.concat_dir} {params.samplename}"

rule plot_CB_coverage:
    """Plot the MT coverage of the single cells"""
    input:
         all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
         barcodes = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_CB_coverage_hist_minReads{num_read}.png"
    shell:
        "python src/plot_CB_coverage.py {input} {output}"


rule scPileup_MT_matrix:
    """Create the position-by-cell coverage matrix"""
    input:
        all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
        barcode_p = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output:
        sc_coverage_f = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
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
        all = "{results}/data/{sample}/MT/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.strands.txt.gz",
        barcode_p = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_barcode_data.p"
    output:
        "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz"
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
        sc_coverage_f = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
    output:
        save_f_heat = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position.png",
        save_f_coverage = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}_MT_position_coverage.png"
    shell:
        #"python src/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_coverage}"
        "python src/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_heat} {output.save_f_coverage}"


def get_filt(w):
    return w.min_cells, w.min_reads, w.topN, w.het_thresh, w.min_het_cells, w.het_count_thresh, w.bq_thresh


rule create_filters:
    input:
        concat_dir = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/{sample}.coverage.strands.txt.gz"
    #output:  "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/af_by_cell.tsv"
    output:
        cov = "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/{sample}.coverage.txt"
    params:
        concat_d = lambda wildcards, input: dirname(input.concat_dir),
        ref_fa = config['mt_ref_fa'],
        name = lambda wildcards: wildcards.sample,
        filt_params = get_filt,
    resources:
        mem_mb=90000
    #log: "{results}/data/{sample}/MT/cellr_{cellr_bc}/{sample}_{num_read}/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/.outcfg"
    shell: "python src/calculate_AF_by_cell.py {params.concat_d} {output.af_f} {params.ref_fa} {params.name} {params.filt_params}"# --log {log}"
