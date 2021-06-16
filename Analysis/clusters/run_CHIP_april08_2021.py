# from src.mlflow.clones import MTExp
# import src.mlflow.mlflow_utils as mu
# #from mlflow.utils import mlflow_tags
# import os
# import mlflow
# import click
# from src.utils.parse_config import read_config_file
from src.vireo.lineage_enrichment import lineage_enrichment


lineage_enrichment("/data2/mito_lineage/data/processed/mttrace/CHIP_april08_2021_Croker/MTblacklist/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/vireoIn/clones",
                   "/data2/mito_lineage/data/processed/mttrace/CHIP_april08_2021_Croker/MTblacklist/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/vireoIn/enrichment",
                   n_clone_list=(20,40), samples="Control,Flt3l,Input", plot_ind=True)


# @click.command()
# @click.argument("entry", type=click.Choice(["main", "preprocess", "results"]))
# @click.argument("configfile", type=click.Path(exists=True))
# @click.option("--preprocmode", default="filters", type=click.Choice(["filters", "vireo"]))
# @click.option("--datauri", default="", type=click.Path())
# def main(configfile, ):
#
#     lineage_enrichment(clones_indir, outdir, n_clone_list, samples,
#                        plot_ind)
#     "python {params.script} {params.donors_indir} {params.clones_indir} {params.OUTDIR} {params.N_DONORS} {params.nclones} {params.samples}"
#
#
# if __name__ == "__main__":
#     main()
    # main(['mtexp_a.yaml',  'main',
    #       '--preproc_mode', 'filters',
    #       '--datauri', '.'])
    # main(['mtexp_a.yaml',  'main',
    #       '--preproc_mode', 'filters',
    #       '--datauri', '.'])
