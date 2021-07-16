#configfile: "parameters/Zhu_Single_Cell_10X_genomics_Human2_002.yaml"
from os.path import join
from src.utils.data_io import sparse_to_cellranger_fragments


configfile: "/data2/mito_lineage/parameters/DUPI_april08_2021/mttrace_mtasnucl.yaml"

cfg_anno = config['annotations']
out = cfg_anno["name"]

rule all:
    input:
        expand("output/data/{out}/{out}.fragments.tsv", out=config['annotations']['name'])
        #f"data/processed/external/{out}/{out}.fragments.tsv"#, out=config['annotations']['name'])

def get_cellr_input(wildcards):
    # peaks_f = f"{indir}/peaks.bed"
    # singlecell_f = f"{indir}/single.bed"
    # mtx_f = ""
    out_d = cfg_anno['datadir']
    return join(out_d, cfg_anno['mtx_f']), join(out_d, cfg_anno['peaks_f']), join(out_d, cfg_anno['cells_f'])

rule peaks_to_fragments:
    input: get_cellr_input
    output: "output/data/{out}/{out}.fragments.tsv"
    run: sparse_to_cellranger_fragments(input[0], input[1], input[2], out_f=output[0], n_cpus=24)

rule project_annotation:
    input: "/data2/mito_lineage/data/processed/mttrace/{prefix}/MTblacklist/{sample}/MT/cellr_True/{sample}_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/"
    output:  "/data2/mito_lineage/Analysis/annotation/output/data/{sample}"
    params:
        name= lambda wildcards: wildcards.sample,
        script= join("src", "R_scripts", "annotations", "annotations.ipynb")
    shell: "papermill -p name {params.name} -p params.indir {params.indir} -p outdir {params.outdir} -in_mgatk {params.in_mgatk} {output[0]}"
