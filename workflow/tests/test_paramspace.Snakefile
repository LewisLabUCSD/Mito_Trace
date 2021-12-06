import pandas as pd
from src.utils.parse_config import read_config_file
from os.path import dirname, join
from src.config import ROOT_DIR
from snakemake.utils import Paramspace
#prefix = config["prefix"]
global_cfg = read_config_file(config['config'])

params = pd.read_csv("parameters/pipeline/params_test.csv").astype("object")#: "int"})

#(dataframe, filename_params=None, param_sep='~')

pspace_mt = Paramspace(params.loc[:,global_cfg['mtpreproc']["params"].keys()],
                       filename_params=list(global_cfg['mtpreproc']["params"].keys()))

samples = ["A","B"]

pspace_samples=Paramspace(pd.DataFrame(["A","B"], columns=["samples"]))
print(pspace_samples)
pspace_filt = Paramspace(params.loc[:,global_cfg['filters']["params"].keys()],
                       filename_params=list(global_cfg['filters']["params"].keys()))
pspace_mgatk = Paramspace(params.loc[:,global_cfg['mgatk']["params"].keys()],
                       filename_params=list(global_cfg['mgatk']["params"].keys()))
#pspace_mult = Paramspace(params, params_ft[ filename_params=global_cfg['multiplex']["params"].keys())
pspace_cln = Paramspace(params.loc[:,["method"]],
                       filename_params=["method"])

pspace_var = Paramspace(params.loc[:,["variants"]],
                       filename_params=["variants"])

#Paramspace(params, params_ft[ filename_params=global_cfg['clones']["simple"]["params"].keys())

rule all:
    input:
        expand("output/tests/pspace/{pspace_mt}/{sample}/mt.txt",
               sample=pspace_samples.instance_patterns,
               pspace_mt=pspace_mt.instance_patterns),
        expand("output/tests/pspace/{pspace_mt}/{sample}/{pspace_filt}/filt.txt",
                sample=pspace_samples.instance_patterns,
                pspace_mt=pspace_mt.instance_patterns,
                pspace_filt=pspace_filt.instance_patterns),
        expand("output/tests/pspace/{pspace_mt}/merged/{pspace_filt}/{pspace_cln}/clone.txt",
                pspace_mt=pspace_mt.instance_patterns,
                pspace_filt=pspace_filt.instance_patterns,
                pspace_cln=pspace_cln.instance_patterns)

rule a:
    output: f"output/tests/pspace/{pspace_mt.wildcard_pattern}/{pspace_samples.wildcard_pattern}/mt.txt",
    shell: "touch {output}"


rule b:
    input:
        f"output/tests/pspace/{pspace_mt.wildcard_pattern}/{pspace_samples.wildcard_pattern}/mt.txt"
    output: f"output/tests/pspace/{pspace_mt.wildcard_pattern}/{pspace_samples.wildcard_pattern}/{pspace_filt.wildcard_pattern}/filt.txt",
    params:
        pspace_filt.instance
    shell: "echo {params[0][mincells]} && echo {params[0][hetthresh]} &&  touch {output} "
    # run:
    #     print('filt params')
    #     print(params)
    #     print('hetthresh')
    #     print(params[0]['hetthresh'])
    #     shell("touch {output}")


rule c:
    input:
        expand("output/tests/pspace/{pspace_mt}/{pspace_samples}/{pspace_filt}/filt.txt",
              pspace_samples=pspace_samples.instance_patterns, pspace_mt=pspace_mt.wildcard_pattern,
               pspace_filt=pspace_filt.wildcard_pattern)
        # expand("output/tests/pspace/{{pspace_mt}}/{pspace_samples}/{{pspace_filt}}/filt.txt",
        #           pspace_samples=pspace_samples.instance_patterns)

    output:
        f"output/tests/pspace/{pspace_mt.wildcard_pattern}/merged/{pspace_filt.wildcard_pattern}/{pspace_cln.wildcard_pattern}/clone.txt"
    shell: "touch {output}"
