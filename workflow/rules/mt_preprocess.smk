import pandas as pd
import pickle
from os.path import basename, dirname
from snakemake.utils import min_version
min_version("6.0")
import subprocess as subp
from os.path import join, dirname
from src.config import ROOT_DIR

# wildcard_constraints:
#     output = '((?!coverage_merged))'
#mtpreproc_extractCB_from_bam_03

res = config['results']
samples = config["samples"] #pd.read_table(config["samples_meta"], dtype=str,sep=',').set_index(["sample_name"], drop=False)

# print('samples')
# print(samples)
mt = config['mtpreproc']
mt_parms = mt["params"]
ft = config['filters']
ft_parms = ft["params"]



######################################################################
## Extract the MT and process for the relevant reads with the "CB" cell barcodes
######################################################################
def get_sample_bam(wildcards):
    bam = samples.loc[wildcards.sample,"bam_f"]
    #print(bam)
    return bam


def get_sample_barcodes(wildcards):
    return samples.loc[wildcards.sample, "barcode_f"]


# rule link_bam:
#     input: get_sample_bam
#     output: "{output}/data/{sample}/00_bam/{sample}.bam"
#     shell: 'ln -sr {input} {output}'

rule link_bam:
    #input: get_sample_bam
    input: get_sample_bam
    output: temp("{output}/data/{sample}/00_bam/{sample}.bam")
    shell: 'ln -sr {input} {output}'

rule index_bam:
    """Index the bam file"""
    input: ("{output}/data/{sample}/00_bam/{sample}.bam") #rules.link_bam.output
    output: temp("{output}/data/{sample}/00_bam/{sample}.bam.bai")
    shell: "samtools index {input}"

rule MT_map:
    """Extract the MT genome"""
    input:
        bam = "{output}/data/{sample}/00_bam/{sample}.bam",
        bai = ("{output}/data/{sample}/00_bam/{sample}.bam.bai"),
    output:
        mt_bam = temp("{output}/data/{sample}/MT/{sample}.MT.bam"),
        mt_bai = temp("{output}/data/{sample}/MT/{sample}.MT.bam.bai")
    params:
        bam = lambda wildcards, input: input.bai.replace('.bai', ''),
        mt_chr= mt["mito_character"],

    run:
        shell("samtools view -b {params.bam} {params.mt_chr} > {output.mt_bam}")
        shell("samtools index {output.mt_bam}")


rule get_bigwig:
    """Extract the MT genome"""
    input: "{output}/data/{sample}/MT/{sample}.MT.bam"
    output:
        coverage="{output}/data/{sample}/MT/{sample}.MT.bw"
    shell: "bamCoverage -b {input} -o {output}"


rule barcode_data:
    """Loop through the bam file and extract the barcode information."""
    input:
        mt_bam="{output}/data/{sample}/MT/{sample}.MT.bam",
        mt_bai="{output}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: "{output}/data/{sample}/MT/{sample}_barcode_data.p",
    params: mt_chr=mt["mito_character"]
    shell:
         "python src/mtpreproc/bam_barcodes_function.py {input.mt_bam} {output} {params.mt_chr}"


rule barcode_filter:
    input:
        barcode_p = "{output}/data/{sample}/MT/{sample}_barcode_data.p",
        cellr_f = get_sample_barcodes
    output:
        "{output}/data/{sample}/MT/cellr_{cellrbc}/{sample}_barcode_data.p",
        "{output}/data/{sample}/MT/cellr_{cellrbc}/{sample}_barcode_data.txt"
    params: cellrbc =  True #lambda wildcards: wildcards.cellrbc
    shell: "python src/mtpreproc/filter_barcodes.py {input} {params} {output[0]}"


######################################################################
## Create cell specfic files and then merge them into a pileup matrix
######################################################################
rule extractCB_from_bam_03:
    input:
        "{output}/data/{sample}/MT/{sample}.MT.bam",
        cells="{output}/data/{sample}/MT/cellr_True/{sample}_barcode_data.txt"
    output:
        sam = temp("{output}/data/{sample}/MT/filtered.sam"),
    params:
        cmd = lambda wildcards, input: f"{{ samtools view -H {input[0]} & samtools view {input[0]} | LC_ALL=C grep -F -f {input.cells}; }}"
    run: shell("{params.cmd} > {output.sam} 2> {output.sam}.log")
    #shell:  "{{ samtools view -H {input[0]} & samtools view {input[0]} | LC_ALL=C grep -F -f {input.cells}; }} > {output.sam} 2> {output.sam}.log"

# rule extractCB_from_bam_01:
#     input:
#         "{output}/data/{sample}/MT/{sample}.MT.bam",
#     output:
#         head=temp("{output}/data/{sample}/MT/SAM_HEADER"),
#     shell: "samtools view -H {input[0]} > {output.head}"

# rule extractCB_from_bam_02:
#     """Extract only the barcodes given"""
#     input:
#         "{output}/data/{sample}/MT/{sample}.MT.bam",
#         cells="{output}/data/{sample}/MT/cellr_True/{sample}_barcode_data.txt"
#     output:
#         body = temp("{output}/data/{sample}/MT/filtered_SAM_body"),
#     shell: "samtools view {input[0]} | LC_ALL=C grep -F -f {input.cells} > {output.body}"
#
# rule extractCB_from_bam_03:
#     input:
#         head_f = "{output}/data/{sample}/MT/SAM_HEADER",
#         body = "{output}/data/{sample}/MT/filtered_SAM_body"
#     output:
#         sam = temp("{output}/data/{sample}/MT/filtered.sam")
#     shell: "samtools view -b {input.body} > {output.sam}"

rule extractCB_from_bam_04:
    input:
        sam = "{output}/data/{sample}/MT/filtered.sam"
    output:
        bam = temp("{output}/data/{sample}/MT/filtered.bam"),
    shell: "samtools view -b {input.sam} > {output.bam}"

    
# rule extractCB_from_bam:
#     input:
#         "{output}/data/{sample}/MT/{sample}.MT.bam",
#         "{output}/data/{sample}/MT/{sample}.MT.bam.bai",
#         cells="{output}/data/{sample}/MT/cellr_True/{sample}_barcode_data.txt"
#     output:
#         head = temp("{output}/data/{sample}/MT/SAM_HEADER"),
#         body = temp("{output}/data/{sample}/MT/filtered_SAM_body"),
#         sam = temp("{output}/data/{sample}/MT/filtered.sam"),
#         bam = temp("{output}/data/{sample}/MT/filtered.bam"),
#     # params:
#     #     config["samtools_path"]
#     run:
#          shell("export BAM_FILE='{input[0]}'"),
#          shell("samtools view -H {input[0]} > {output.head}"),
#          shell("samtools view {input[0]} | LC_ALL=C grep -F -f {input.cells} > {output.body}"),
#          shell("cat {output.head} {output.body} > {output.sam}"),
#          shell("samtools view -b {output.sam} > {output.bam}"),
#          shell("samtools index {output.bam}")

rule sortCB:
    input:
        mt_bam="{output}/data/{sample}/MT/filtered.bam",
    output: temp("{output}/data/{sample}/MT/{sample}.MT.CB.bam")
    shell: "samtools sort -t CB {input.mt_bam} > {output}"



rule scBam:
    """Extract each single-cell and put into respective bam file"""
    input: "{output}/data/{sample}/MT/{sample}.MT.CB.bam"
        #mt_bam="{output}/data/{sample}/MT/{sample}.MT.bam",
        #mt_bai="{output}/data/{sample}/MT/{sample}.MT.bam.bai"
    output: temp(directory("{output}/data/{sample}/MT/{sample}_scBam"))
    threads: 18
    #directory("/data/isshamie/miito_lineage/mttrace/{sample}/MT/{sample}_scBam")
    shell:
        "python src/mtpreproc/split_by_CB.py {input} {output}"


rule scPileup:
    """Run the first part of the MT-genotype function by getting the read pileups for each bam file for each nucleotide and overall coverage"""
    input:
        scBam = "{output}/data/{sample}/MT/{sample}_scBam",
        #scBam = "{output}/data/{sample}/MT/{sample}_scBam", "{output}/data/{sample}/MT/cellr_{cellrbc}/{sample}_barcode_data.txt"
        barcodes = "{output}/data/{sample}/MT/{sample}_barcode_data.p",
    output:
          temp(directory("{output}/data/{sample}/MT/{sample}_scPileup_{num_read}"))
    params:
        base_quality = mt_parms['basequality'], # mean base quality threshold for each position-cell, otherwise 0.
        num_reads = lambda wildcards: wildcards.num_read,
        #align_quality = lambda wildcards: wildcards.alignQual
    threads: 18
    log: "{output}/logs/{sample}/MT/{sample}_scPileup_{num_read}/log.txt"
    shell:
         "python src/mtpreproc/scPileup_counts.py {input.scBam} {output} {input.barcodes} {params.num_reads} {params.base_quality} > {log}"


def concat_files(directory, samplename, nt):
    cmd = f"find {directory} -type f -name *.{nt}.txt -exec cat {{}} \\; > {samplename}_all.{nt}.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)
    cmd = f"find {directory} -type f -name *.{nt}.minus.txt -exec cat {{}} \\; > {samplename}_all.{nt}.minus.txt"
    print(cmd)
    subp.check_call(str(cmd), shell=True)

    # # Gzip the files
    # cmd = f"gzip {samplename}_all.{nt}.txt"
    # print(cmd)
    # subp.check_call(str(cmd), shell=True)
    # cmd = f"gzip {samplename}_all.{nt}.minus.txt"
    # print(cmd)
    # subp.check_call(str(cmd), shell=True)
    return


rule scPileup_concat:
    """ Concat the pileups for each separate cell pileup file
    """
    input:
        scPileup_dir = "{output}/data/{sample}/MT/{sample}_scPileup_{num_read}",
        #scBam = "{output}/data/{sample}/MT/{sample}_scBam"
    output:
        temp(expand("{{output}}/data/{{sample}}/MT/scPileup_concat_{{num_read}}/numread_{{num_read}}_all.{nt}{strand}.txt",
             nt=["coverage", "A", "C", "G", "T"], strand=[".minus", ""]))
    threads: 28
    params:
        samplename = lambda wildcards, output: output[0].split("_all.coverage.minus.txt")[0]
          #"data/processed/{sample}/scPileup_concat/numread_{num_read}"
    run:
        for n in ["A", "C", "G", "T", "coverage"]:
            concat_files(input.scPileup_dir, params.samplename, n)


rule scPileup_mergeStrands:
    """ Run the second part of the MT-genotype pipeline, which just concatenates all the pileup data for each nucleotide and overall."""
    input:
        expand("{{output}}/data/{{sample}}/MT/scPileup_concat_{{num_read}}/numread_{{num_read}}_all.{nt}{strand}.txt",
               nt=["coverage", "A", "C", "G", "T"], strand=[".minus", ""])
    output:
        all = (expand("{{output}}/data/{{sample}}/MT/scPileup_concat_{{num_read}}/numread_{{num_read}}_all.{nt}.strands.txt",
                     nt=["coverage", "A", "C", "G", "T"]))
    params:
        concat_dir = lambda wildcards, input: dirname(input[0]),
        samplename = lambda wildcards, output: basename(output.all[0]).split("_all.coverage.strands.txt")[0]
          #"data/processed/{sample}/scPileup_concat/numread_{num_read}"
    log: "{output}/logs/{sample}/MT/scPileup_concat_{num_read}/log.txt"
    shell:
         "python src/mtpreproc/scPileup_mergeStrands.py {params.concat_dir} {params.samplename} > {log}"


rule filter_cell_bc:
    """Extracts only the relevant cell barcodes and removes the .bam from the barcode names."""
    input:
        all = expand("{{output}}/data/{{sample}}/MT/scPileup_concat_{{num_read}}/numread_{{num_read}}_all.{nt}.strands.txt",
                     nt=["coverage", "A", "C", "G", "T"]),
        barcode_p = "{output}/data/{sample}/MT/cellr_{cellrbc}/{sample}_barcode_data.p"
    output:
        ("{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt")
    run:
        for n in ["A", "C", "G", "T", "coverage"]:
            #print('n', n)
            curr_f = input.all[0]
            curr_f = curr_f.replace(".coverage.", "." + n + ".")
            #print('curr_f', curr_f)
            df = pd.read_csv(curr_f, sep=',')
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
            df.to_csv(curr_out_f, header=None, index=None)


rule scPileup_MT_matrix:
    """Create the position-by-cell coverage matrix"""
    input:
        all = rules.scPileup_mergeStrands.output.all[0],  #"{output}/data/{sample}/MT/scPileup_concat_{num_read}/numread_{num_read}_all.coverage.strands.txt.gz",
        barcode_p = rules.barcode_filter.output[0] #"{output}/data/{sample}/MT/cellr_{cellrbc}/{sample}_barcode_data.p"
    output:
        sc_coverage_f = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/sc_coverage.csv"
    threads: 16
    params:
        maxBP= 16571
    shell:
        "python src/figures/plot_heatmap_coverage.py sc_mt {input.barcode_p} {input.all} {output.sc_coverage_f} {params.maxBP}"


rule plot_sc_coverageBar_and_heat:
    """Plot the posiitonal coverages."""
    input:
        sc_coverage_f = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/sc_coverage.csv"
    output:
        save_f_coverage = report("{output}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position_coverage.png", category="mtpreproc"),
        save_f_heat = report("{output}/figures/{sample}/MT/cellr_{cellrbc}/numread_{num_read}_MT_position.png", category="mtpreproc"),
    shell:
        "python src/figures/plot_heatmap_coverage.py plot {input.sc_coverage_f} {output.save_f_heat} {output.save_f_coverage}"


def get_filt(w):
    return w.mincells, w.minreads, w.topN, w.hetthresh, w.minhetcells, w.hetcountthresh, w.bqthresh




# rule create_filters:
#     input:
#         concat_dir = (rules.filter_cell_bc.output[0]) #"{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt.gz"
#     output:
#         cov = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/{sample}.coverage.txt",
#         af = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/af_by_cell.tsv",
#         fig = report("{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/heatmap.png",
#                       category="mtpreproc", subcategory="filter"),
#         fig2 = report("{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/initial_cell_depth.png",
#                       category="mtpreproc",subcategory="filter"),
#     params:
#         concat_d = lambda wildcards, input: dirname(input.concat_dir),
#         ref_fa = config["genome_path"][config['genome']]['mt_ref_fa'],
#         name = lambda wildcards: wildcards.sample,
#         filt_params = get_filt,
#     resources:
#         mem_mb=90000
#     shell: "python src/mtpreproc/calculate_AF_by_cell.py {params.concat_d} {output.af} {params.ref_fa} {params.name} {params.filt_params}" # --log {log}"
#

rule create_filters_v02:
    input:
        concat_dir = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt" #(rules.filter_cell_bc.output[0]) #"{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/{sample}.coverage.strands.txt.gz"
    output:
        note = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/out.ipynb",
        cov = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/{sample}.coverage.txt",
        af = "{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/af_by_cell.tsv",
        fig = report("{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/heatmap.png",
                      category="mtpreproc", subcategory="filter"),
        fig2 = report("{output}/data/{sample}/MT/cellr_{cellrbc}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/initial_cell_depth.png",
                      category="mtpreproc",subcategory="filter"),
    params:
        concat_d = lambda wildcards, input: dirname(input.concat_dir),
        ref_fa = config["genome_path"][config['genome']]['mt_ref_fa'],
        name = lambda wildcards: wildcards.sample,
        #filt_params = get_filt,
        note = join(ROOT_DIR, "workflow/notebooks/af_filter/filter_af_by_cell_v02.ipynb")
    version: "v02"
    resources:
        mem_mb=90000
    shell:
        """
        papermill -p scpileup_dir {params.concat_d} -p af_f {output.af} -p mt_ref_fa {params.ref_fa} -p name {params.name} \
        -p min_cells {wildcards.mincells} -p min_reads {wildcards.minreads} -p topn {wildcards.topN} \
        -p het_thresh {wildcards.hetthresh} -p min_het_cells {wildcards.minhetcells} -p het_count_thresh {wildcards.hetcountthresh} \
        -p bq_thresh {wildcards.bqthresh} {params.note} {output.note}
        """


