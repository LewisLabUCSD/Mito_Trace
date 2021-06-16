#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"

# Rkernel: /usr/share/R
# rule all:
#     input:
#

rule project_annotation:
    input: "/data2/mito_lineage/data/processed/mttrace/{prefix}/MTblacklist/{sample}/MT/cellr_True/{sample}_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/"
    output:  "/data2/mito_lineage/Analysis/annotation/output/data/{sample}"
    params: name= lambda wildcards: wildcards.sample
    shell: "papermill -p name {params.name} -p params.indir {params.indir} -p outdir {params.outdir} -in_mgatk {params.in_mgatk} {output[0]}"
