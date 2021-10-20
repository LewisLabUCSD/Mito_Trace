from src.utils.loop_smk import wrap_run_snakemake
from src.config import ROOT_DIR
from os.path import join
from src.utils.parse_config import read_config_file
import snakemake


cfg_f= join(ROOT_DIR, "parameters/pipeline/wrappers", "wrap_pipeline.yaml")
#wrap_run_snakemake(cfg_f)

config = read_config_file(cfg_f)  # join(ROOT_DIR,"parameters/pipeline/wrap_annotation.yaml"))
print(config['indir'])
for i in config["indir"]:
    print('i', i)
    curr_config = join(config["params_dir"], config["indir"][i])
    print(curr_config)
    print(config['pipename'])

    out = snakemake.snakemake(config["snakefile"], configfiles=[curr_config],
                              dryrun=False,
                              forcetargets=['knn_enrichment',
                                            'vireo_enrichment'],
                              cores=8, verbose=True)
