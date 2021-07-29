""" Wrapper to run snakemake, create report, and add to git

If some stages are already run, it will pick up from where it left off,
unless overwrite=True
"""
import click
import os
from os.path import basename, dirname, join
import snakemake
from src.utils.parse_config import read_config_file, write_config_file
#from snakemake.utils import para
from src.utils.paramspace import Paramspace
import mlflow
import src.mlflow.mlflow_utils as mu
from icecream import ic
import subprocess


def run_git(paths:list, commit_msg, to_push=False, src="origin", remote="master",
            branch='workflow_results'):
    assert(len(paths)>0)

    # check if branch already made. 0 returncode means it is there
    is_branch = subprocess.run(["git", "show-ref", "--verify", "--quiet", f"refs/heads/{branch}"], universal_newlines=True).returncode
    if is_branch==0:
        os.system(f"git checkout {branch}")
    else:
        os.system(f"git checkout -b {branch}")
    for p in paths:
        cmd = f"git add {p}"
        os.system(cmd)
    cmd = f'git commit -m "{commit_msg}"'
    os.system(cmd)
    if to_push:
        cmd = f'git push -u {src} {branch}'
        os.system(cmd)
    return



# def run(**kwargs):
#     check_required(kwargs, ["protocol"])
#     if kwargs["protocol"] not in PROTOCOLS:
#         raise InvalidParameterError(
#             "Invalid protocol selection: " +
#             "{}. Valid protocols are: {}".format(
#                 kwargs["protocol"], ", ".join(PROTOCOLS.keys())
#             )
#         )
#     return PROTOCOLS[kwargs["protocol"]](**kwargs)
#

def process_params(skeleton):
    return


def setup_files(indir, params, pipe='pipeline'):
    return


def is_mlflow():
    return


def get_curr_dir(params, is_mlflow=True):
    return



def check_input(config, files_config, curr_p, git_commit=None,
                mlflow_tracking_uri=None):
    """ Check if the input for each pipeline exists

    1.
    :param config:
    :param entry_point:
    :param files_config:
    :return:
    """
    if mlflow_tracking_uri is not None:
        mlflow.set_tracking_uri(mlflow_tracking_uri)
    for p in config["stages"]:
        local_indirs = files_config[p].get("input", [])
        external_indirs = files_config[p].get("external", [])

        # See if params were already ran
        existing_run = mu._already_ran(p, config[p]['params'], git_commit)
        if existing_run is None:
            return False

        # If done, add
    return True

@click.command()
@click.argument("smkfile", type=click.Path(exists=True))
@click.argument("configfile", type=click.Path(exists=True))
@click.option("--outdir", default=".", type=click.Path(exists=True))
@click.option("--to_git", default=False)
@click.option("--targets", "-t", type=click.STRING, multiple=True)
@click.option("--dryrun", default=False)
@click.option("--mlflow", default=False)
@click.option("--forcetargets", default=False)
@click.option("--to_gitpush", default=False)
@click.option("--template_cfg", default="")
def main(smkfile, configfile, outdir, to_git, targets, dryrun, mlflow, forcetargets, to_gitpush, template_cfg):
    """ Runs snakemake and/or mlflow pipeline


    OUTPUT Folders:
    {Experiment Name (entrypoint)}/{PREFIX}/{mlflow#}/{folder
    :param smkfile:
    :param configfile:
    :param outdir:
    :param to_git:
    :param target:
    :param dryrun:
    :return:
    """
    ic(targets)
    if len(targets) == 0:
        targets=None
    ic(targets)
    if mlflow:
        mlflow.set_tracking_uri(outdir)

    # Setup config output file names for the target
    config = read_config_file(configfile)
    pipeline = config.get('pipeline', [])
    outdir = join(outdir, config['prefix'])
    ic(configfile)
    ic(outdir)

    # Process each input file to get the full paths
    #env_files = load_env_files(env)

    # Compare parameter file to schema
    #validate_schema_files(parameters, parameters_schema)
    for p in pipeline:
        setup_files(outdir, configfile, pipe=p)
        if "params" in config[p]:
            Paramspace(config[p], filename_params="*", param_sep='_')

    snakemake.snakemake(smkfile, configfiles=[configfile],
                        targets=targets, dryrun=dryrun, forcetargets=forcetargets)

    # Make the report
    report = join(outdir, f"report_{basename(smkfile)}.html" )
    if not dryrun:
        snakemake.snakemake(smkfile, configfiles=[configfile],
                            targets=targets, report=report)
        # copy over snakemake and configfile to output:
        os.system(f"cp {configfile} {outdir}/{basename(configfile)}.incfg")
        os.system(f"cp {smkfile} {outdir}/{basename(smkfile)}.insmk")
        if to_git:
            run_git([outdir], commit_msg=f"Ran pipeline for {basename(smkfile)}, {configfile} saved to {outdir} [results]",
                    to_push=to_gitpush)
        write_config_file(join(outdir, "params.outcfg"),config)


    return


if __name__ == "__main__":
    main()
