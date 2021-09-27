from src.utils.parse_config import read_config_file
from src.utils import run_snakemake


def wrap_run_snakemake(config_f):
    config = read_config_file(config_f)#join(ROOT_DIR,"parameters/pipeline/wrap_annotation.yaml"))
    for i in config["indir"]:
        curr_config = (config["indir"][i])
        run_snakemake.main(smkfile=config['snakefile'],
                           configfile=curr_config,
                           pipename=config['pipename'], outdir=None, to_git=True, targets='all',
                           dryrun=True, mlflow=None,
                           to_gitpush=True, cores=8)

