wildcard_constraints:
    num_read="\d+"

import os
import pandas as pd
from snakemake.utils import validate

#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"

#raw = config["raw"]

#bam = config["bam"]
num_reads_filter = config["num_reads_filter"]
maxBP = config["maxBP"]
ref_fa = config["ref_fa"]
print(pd.read_table(config["samples"], dtype=str,sep=','))
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample"], drop=False)
#(samples)

#print('index',samples.index)
res = config["results"]
mq = config["mapq"]


#workdir: config["work_dir"]


wildcard_constraints:
    mapq="\d+",
    cellr='True|False'

rule all:
    input:
        expand("{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/AF/af_filt_{min_cells}_{min_reads}_{top_pos}_{top_cells}_{cell_mt_coverage}_hetC{het_thresh}_hetC{min_het_cells}.csv", results=res, sample=samples["sample"].values,
               mapq = mq, cellr_bc=config["use_cellr_barcode"],
               num_read=config["num_reads_filter"],
               min_cells=config["min_cells"],
               min_reads=config["min_reads"],
               top_pos = config["top_pos"],
               top_cells = config["top_cells"],
               cell_mt_coverage=config["cell_mt_coverage"],
               het_thresh=config['het_thresh'], min_het_cells=config['min_het_cells']),



#Old af_by_cell (not tested) where
#"""Create the AF-by-cell csv file"""
rule create_AF_by_cell:
    input:
        coverage_dir = "{results}/{sample}/mapq_{mapq}/{sample}_scPileup_{num_read}/",
        concat_dir = "{results}/{sample}/mapq_{mapq}/scPileup_concat_{num_read}/"  #{sample}_{num_read}_all.coverage.txt.gz",
    output:  "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/AF/af_filt_{min_cells}_{min_reads}_{top_pos}_{top_cells}_{cell_mt_coverage}_hetC{het_thresh}_hetC{min_het_cells}.csv"
    params:
        min_cells="{min_cells}",
        min_reads="{min_reads}",
        top_cells ="{top_cells}",
        cell_mt_coverage ="{cell_mt_coverage}",
        het_thresh='{het_thresh}',
        min_het_cells='{min_het_cells}',
        maxBP = maxBP,
        ref_fa = config["ref_fa"]
    shell:
        "python src/calculate_AF_by_cell.py {input} {output} {params}"




#barcode_filter : "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_barcode_data.p"

#
#
#
# """Create the AF-by-cell pickle file"""
# rule generate_allele_frequencies:
#     input:  "{results}/{sample}/mapq_{mapq}/scPileup_concat_{num_read}/{sample}_{num_read}_all.coverage.txt.gz"
#     output:  "{results}/{sample}/mapq_{mapq}/AF_min{num_read}.tensor.p"
#     params: maxBP = config['maxBP']
#     shell:
#       "python src/calculate_AF.py {input} {output} -maxBP {params}"
#
#
#
#
#
# # position_filter: # At least n cells with greater than n reads
# #   min_cells:
# #     - 100
# #     - 500
# #     - 10
# #   min_reads:
# #     - 100
# #
# # cell_mt_coverage:
# #   - 1
# #   - 10
# #   - 100
#
#         # min_cells = ,
#         # min_reads = ,
#         # cell_mt_coverage=
# #position_bq_thresh:
# """Create a text file of cell barcodes to keep"""
# rule filter_af:
#     input:
#          af_f = "{results}/{sample}/mapq_{mapq}/AF_min{num_read}.tensor.p",
#          barcode_data = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_barcode_data.p"
#         #"{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/sc_coverage.csv"
#     output:
#         af_filt_f = "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/af_filt_{min_cells}_{min_reads}_{top_pos}_{top_cells}_{cell_mt_coverage}.csv"
#     params:
#         min_cells = "{min_cells}",
#         min_reads = "{min_reads}",
#         top_cells = "{top_cells}",
#         top_pos = "{top_pos}",
#         cell_mt_coverage = "{cell_mt_coverage}",
#
#     shell: "python src/af_filters.py {input} {output} {params}"
#
#
# rule aggregate_samples:
#     input:
#         expand("{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/af_filt_{min_cells}_{min_reads}_{top_pos}_{top_cells}_{cell_mt_coverage}.csv", results=res, sample=samples["sample"].values,
#                mapq = mq, cellr_bc=config["use_cellr_barcode"],
#                num_read=config["num_reads_filter"],
#                min_cells=config["min_cells"],
#                min_reads=config["min_reads"],
#                top_pos = config["top_pos"],
#                top_cells = config["top_cells"],
#                cell_mt_coverage=config["cell_mt_coverage"])
#     output: "{results}/mq_{mapq}_cellr_{cellr_bc}_nr_{num_reads}_mc_{min_cells}_mr_{min_reads}/aggregate_af.csv"
#     shell: "python aggregate.py {input} {output}"
#
# rule plot_cells:
#     input: "{results}/mq_{mapq}_cellr_{cellr_bc}_nr_{num_reads}_mc_{min_cells}_mr_{min_reads}/aggregate_af.csv"
#     output: "{results}/{sample}/mapq_{mapq}/cellr_{cellr_bc}/{sample}_{num_read}/cluster_{min_cells}_{min_reads}_{cell_mt_coverage}.png"
#     shell: "python plot_lineage.py {input} {output}"

#
# rule filter_positions:
#""" Create a text file of variants to keep
#     input:
#     output:
#     shell:
#
# rule generate_allele_frequencies:
#"""Create the AF-by-cell csv file"""
#     input:
#     output:
#     shell:


#rule plot_allele_frequencies:
#"""Plot the AF-by-cell heatmap"""
