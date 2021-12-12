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
import src.mlflow.mlflow_utils as mu
import mlflow
from icecream import ic
import subprocess
import datetime
import pandas as pd


def run_git(paths:list, commit_msg, to_push=False, src="origin", remote="master",
            branch='master'):
    assert(len(paths)>0)

    # check if branch already made. 0 returncode means it is there
    is_branch = subprocess.run(["git", "show-ref", "--verify", "--quiet", f"refs/heads/{branch}"], universal_newlines=True).returncode
    try:
        if is_branch==0:
            try:
                os.system(f"git checkout {branch}")
            except:
                print(f"Already on branch {branch}")
        else:
            os.system(f"git checkout -b {branch}")
    except:
        print("Failed to checkout git branch")
        return
    try:
        for p in paths:
            cmd = f"git add {p}"
            os.system(cmd)
    except:
        print("Failed to add to git")
        return

    cmd = f'git commit -m "{commit_msg}"'
    os.system(cmd)

    if to_push:
        try:
            cmd = f'git push -u {src} {branch}'
            os.system(cmd)
        except:
            print("Failed to push to git")
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


def create_params_combination(params, subset=None):
    """ Creates all combinations of the different parameters given and puts
    into long pandas df. This can then be input to snakemake

    :param params: [{param:[param instances]}] list of the parameters to use and the instances
    :param subset: subset of the params keys to use
    :return:
    """
    from itertools import product
    if subset is not None:
        params = [list(p.keys())[0] for p in params if list(p.keys())[0] in subset]

    params_df = pd.DataFrame(product([i.values() for i in params]),
                    columns=[list(i.keys())[0] for i in params])
    return params_df


def create_params_combination_global(cfg):
    params = {}
    if 'pipeline' in cfg:
        analyses = cfg["pipeline"]
        for an in analyses:
            params[an] = analyses["params"]
    params_df = create_params_combination(params)
    return params_df


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


def run(smkfile, configfile, pipename, outdir=None, to_git=False, targets=None, dryrun=False, forcetargets=False,
        to_gitpush=False, cores=2):
    ic(targets)
    if targets==None or len(targets) == 0:
        targets=None
    ic(targets)

    curr_time = datetime.datetime.now().strftime("%B_%d_%Y_%H%M%S")

    out = snakemake.snakemake(smkfile, configfiles=[configfile],
                        targets=targets, dryrun=dryrun,
                              forcetargets=forcetargets,
                              cores=cores, force_incomplete=True)
    print("out")
    print(out)

    # Setup config output file names
    config = read_config_file(configfile)
    #pipeline = config.get('pipeline', [])
    if 'outdir' in config:
        #config['outdir']
        if not os.path.exists(join(config["outdir"], pipename)):
            raise OSError(f'Expected dir {join(config["outdir"], pipename)} but does not exist')
        outdir = join(config['outdir'], pipename, config['prefix'], "configs")
    else:
        if not os.path.exists(join(outdir, pipename)):
            raise OSError(f'Expected dir {join(outdir, pipename)} but does not exist')
        outdir = join(outdir, pipename, config['prefix'], "configs")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ic(configfile)
    ic(outdir)


    if not dryrun:
        # Make the report
        report_f = join(outdir,
                      f"report_{basename(smkfile)}_cfg_{basename(configfile)}_{curr_time}.html")
        ic(report_f)
        snakemake.snakemake(smkfile, configfiles=[configfile],
                            targets=targets, report=report_f,
                            force_incomplete=True)

        ################
        # copy over snakemake and configfile to output:
        ################
        incfg_f = f"{outdir}/{basename(configfile)}_{curr_time}.incfg"
        insmk_f = f"{outdir}/{basename(configfile)}_{curr_time}.insmk"
        params_out_f = join(outdir, f"params.outcfg _{curr_time}")
        os.system(f"cp {configfile} {incfg_f}")
        os.system(f"cp {smkfile} {insmk_f}")
        write_config_file(params_out_f, config)

        snakemake.snakemake(smkfile, configfiles=[configfile],
                            targets=targets,
                            printrulegraph=True)
        cmd = f"snakemake -s {smkfile} --configfile {configfile} --dag | dot -Tsvg > {outdir}/dag_{curr_time}.svg"
        print(cmd)
        os.system(cmd)
        # if to_git:
        #     ## TODO
        #     run_git([report_f,
        #              params_out_f,
        #              insmk_f,
        #              incfg_f],
        #             commit_msg=f"Ran pipeline for {basename(smkfile)}, {configfile} saved to {outdir} [results]",
        #             to_push=to_gitpush)


        # 3. Update _stateOfAnalysis.txt file
        an_f = join(outdir, '_stateOfAnalysis.txt')
        # Timestamp and update
        s = 'Ran and Needs inspection\n' + 'Time:' + curr_time + "\n"
        with open(an_f, 'a') as f:
            f.write(s)
    return


@click.command()
@click.argument("smkfile", type=click.Path(exists=True))
@click.argument("configfile", type=click.Path(exists=True))
@click.option("--pipename", "-p", default="pipeline", type=click.STRING)
@click.option("--outdir", "-o", default=".", type=click.Path(exists=True))
@click.option("--to_git", default=False)
@click.option("--targets", "-t", type=click.STRING, multiple=True)
@click.option("--dryrun", default=False)
@click.option("--forcetargets", default=False)
@click.option("--to_gitpush", default=False)
@click.option("--cores", default=4)
def main(smkfile, configfile, pipename, outdir, to_git, targets, dryrun, forcetargets, to_gitpush, cores):
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

    run(smkfile, configfile, pipename, outdir, to_git, targets, dryrun,
        forcetargets, to_gitpush, cores)

    return


if __name__ == "__main__":
    main()
