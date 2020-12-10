#import pymc3 as pm
import numpy as np
from numpy import random
import os
import pandas as pd
import pickle
from src.simulations.utils.config import read_config_file, write_config_file
from src.simulations.utils.config import check_required

import click

@click.group(help="Pick a command")
@click.argument('f_save', type=click.Path())
@click.argument('prefix', type=click.STRING)
@click.pass_context
def main(ctx, f_save, prefix):
    """
    Read in parameters and set defaults
    :return:
    """
    ctx.ensure_object(dict)
    ctx.obj['f_save'] = f_save
    ctx.obj['prefix'] = prefix
    return


@main.command()
@click.argument('p', default='hi')
@click.option('--protocol', default='protocol')
@click.pass_context
def run_one_sim(ctx, p, protocol):
    f_save = ctx.obj['f_save']
    # if protocol == "main":
    #     check_required(p, ['f_save', 'prefix'])
    # elif protocol == "binomial":
    #     check_required(p['binomial'], ['steps'])
    print('f_save', f_save)
    print(p)
    print(protocol)

@main.command()
@click.argument('sweep_params', type=click.Path())
@click.option('--steps', default=10)
@click.pass_context
def run_hyperparams(ctx, sweep_params, steps):
    f_save = ctx.obj['f_save']
    print('f_save', f_save)
    # if protocol == "main":
    #     check_required(p, ['f_save', 'prefix'])
    # elif protocol == "binomial":
    #     check_required(p['binomial'], ['steps'])
    print('steps', steps)
    print("sweep_params", sweep_params)

def grow_parameters():
    return


if __name__ == "__main__":
    main(obj={})

    # main(['f_save', "foo.txt",
    #       "prefix", "foo"])
