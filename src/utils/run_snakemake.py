""" Wrapper to run snakemake, create report, and add to git

If some stages are already run, it will pick up from where it left off,
unless overwrite=True
"""
import click
import os
import snakemake
from src.utils.parse_config import read_config_file
#from snakemake.utils import para
from src.utils.paramspace import Paramspace

def run_git(paths:list, commit_msg, to_push=False, src="origin", remote="master"):
    assert(len(paths)>0)
    for p in paths:
        cmd = f"git add {p}"
        os.system(cmd)
    cmd = f'git commit -m "{commit_msg}"'
    os.system(cmd)
    if to_push:
        cmd = f'git push {src} {remote}'
        os.system(cmd)
    return



@click.argument("smkfile", type=click.Path(exists=True))
@click.argument("configfile", type=click.Path(exists=True))
@click.argument("outdir", type=click.Path(exists=True))
@click.argument("to_git", type=click.BOOL)
@click.option("--target", default=None, type=click.STRING)
@click.option("--dryrun", default=True, type=click.STRING)
def main(smkfile, configfile, outdir, to_git, target, dryrun):
    # Setup config output file names for the target
    config = read_config_file(configfile)
    snakemake.snakemake(smkfile, configfiles=[configfile],
                        targets=target, dryrun=dryrun)
    report = os.path.join(outdir, "report.html" )
    for p in config['pipeline']:
        Paramspace(config[p], filename_params="*", param_sep='_')
    snakemake.snakemake(smkfile, configfiles=[configfile],
                        targets=target, report=report)
    # copy over snakemake and configfile to output:
    os.system(f"cp {configfile} {outdir}/{os.path.dirname(configfile)}.incfg")
    os.system(f"cp {smkfile} {outdir}/{os.path.dirname(smkfile)}.insmk")
    if to_git:
        run_git([outdir], commit_msg=f"Ran pipeline for {smkfile}, {configfile} saved to {outdir}")
    return

if __name__ == "__main__":
    main()
