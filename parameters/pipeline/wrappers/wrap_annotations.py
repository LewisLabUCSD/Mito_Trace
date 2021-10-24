from src.utils.loop_smk import wrap_run_snakemake
from src.config import ROOT_DIR
from os.path import join
from src.utils.parse_config import read_config_file
import snakemake


cfg_f= join(ROOT_DIR, "parameters/pipeline/wrappers", "wrap_annotations.yaml")
wrap_run_snakemake(cfg_f, dryrun=False)

# config = read_config_file(cfg_f)  # join(ROOT_DIR,"parameters/pipeline/wrap_annotation.yaml"))
# print(config['indir'])
# for i in config["indir"]:
#     print('i', i)
#     curr_config = join(config["params_dir"], config["indir"][i])
#     print(curr_config)
#     print(config['pipename'])
#
#    # configfiles = [join(config["params_dir"], x) for x in config["indir"].values()]
#     #print('configfiles', configfiles)
#     out = snakemake.snakemake(config["snakefile"], configfiles=[curr_config],
#                               dryrun=False,
#                               forcetargets='createMergedExpSignac',
#                               cores=8, verbose=True)
