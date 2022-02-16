from src.utils.loop_smk import wrap_run_snakemake
from src.config import ROOT_DIR
from os.path import join
from src.utils.parse_config import read_config_file
import snakemake
from src.utils.run_snakemake import run


cfg_f = join(ROOT_DIR, "parameters/pipeline/wrappers", "wrap_pipeline.yaml")
#wrap_run_snakemake(cfg_f)

config = read_config_file(cfg_f)  # join(ROOT_DIR,"parameters/pipeline/wrap_annotation.yaml"))
print(config['indir'])
if "targets" not in config:
    targets = None
    forcetargets=False
else:
    targets = config["targets"]
    forcetargets = True
    print('targets', targets)

dryrun = config.get("dryrun", False)
to_report = config.get("to_report", False)

for i in config["indir"]:
    print('i', i)
    curr_config = join(config["params_dir"], config["indir"][i])
    print(curr_config)
    print(config['pipename'])

    # out = snakemake.snakemake(config["snakefile"], configfiles=[curr_config],
    #                           dryrun=False,
    #                           forcetargets=['knn_enrichment',
    #                                         'vireo_enrichment'],
    #                           cores=8, verbose=True)
    run(config["snakefile"], curr_config, config["pipename"], outdir=None,
        to_git=False, targets=None, dryrun=dryrun, forcerun=targets, to_report=to_report,
            forcetargets=forcetargets, to_gitpush=False, cores=16)