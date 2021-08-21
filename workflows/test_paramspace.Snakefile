from snakemake.utils import Paramspace
import pandas as pd
from src.utils.parse_config import read_config_file
from os.path import dirname, join
from src.config import ROOT_DIR
prefix = config["prefix"]
global_cfg = read_config_file(config['config'])

params = pd.read_csv(config["params"]) #parameters/pipeline/pipeline.yaml


pspace_mt = Paramspace(params, params_ft[ filename_params=global_cfg['mtpreproc']["params"].keys())
pspace_filt = Paramspace(params, params_ft[ filename_params=global_cfg['filters']["params"].keys())
pspace_mgatk = Paramspace(params, params_ft[ filename_params=global_cfg['mgatk']["params"].keys())
#pspace_mult = Paramspace(params, params_ft[ filename_params=global_cfg['multiplex']["params"].keys())
pspace_cln = Paramspace(params, params_ft[ filename_params=global_cfg['clones']["simple"]["params"].keys())


# rule all:
#     input:
#         expand("output/tests/pipeline/{prefix}/results/merged/{pspace_mt}/{pspace_filt}/dendrograms/af_dendro.ipynb",
#                prefix=prefix, pspace_mt=pspace_mt.instance_patterns, pspace_filt=pspace_filt.instance_patterns),
#         # expand("output/tests/pipeline/{prefix}/figures/merged/{mt_p}/{filter_p}/{mgatk_p}/clones/{clones_params}/variants.ipynb",
#         #         prefix=config["prefix"],
#         #         mt_p = pspace_mt.instance_patterns,
#         #         filter_p = pspace_filt.instance_patterns,
#         #         mgatk_p = pspace_mgatk.instance_patterns,
#         #         clones_params=pspace_cln.instance_patterns) #config['multiplex']['n_clone_list'],
#
# rule a:
#     output:
#          f"output/tests/pipeline/{prefix}/results/merged/{pspace_mt.wildcard_pattern}/{pspace_filt.wildcard_pattern}/dendrograms/af_dendro.ipynb"
#     params:
#         mt = pspace_mt.instance #pspace_cln.instance.get("n_clones")
#     shell: "echo {params.mt} && touch {output}"


# rule b:
#     output:
#         expand(f"output/tests/pipeline/{prefix}/figures/merged/{pspace_mt.wildcard_pattern}/{pspace_filt.wildcard_pattern}/{pspace_mgatk.wildcard_pattern}/clones/{pspace_cln.wildcard_pattern}/variants.ipynb"),
#     #params:
#     shell:
#          "touch {output}"


# rule all:
#     input:
#         expand("output/tests/pipeline/{prefix}/results/merged/{pspace_mt}/dendrograms/af_dendro.ipynb",
#                prefix=prefix, pspace_mt=pspace_mt.instance_patterns)
        # expand("output/tests/pipeline/{prefix}/figures/merged/{mt_p}/{filter_p}/{mgatk_p}/clones/{clones_params}/variants.ipynb",
        #         prefix=config["prefix"],
        #         mt_p = pspace_mt.instance_patterns,
        #         filter_p = pspace_filt.instance_patterns,
        #         mgatk_p = pspace_mgatk.instance_patterns,
        #         clones_params=pspace_cln.instance_patterns) #config['multiplex']['n_clone_list'],

rule all:
    input:
        expand("output/tests/pipeline/{prefix}/results/{sample}/output.txt",
               prefix=prefix, sample=['A','B'])
rule a:
    input:
        "output/tests/pipeline/{prefix}/foo.txt",
        cov =  "output/tests/pipeline/{prefix}/bar.txt",
    output:
         expand("output/tests/pipeline/{{prefix}}/results/{sample}/output.txt",
                 prefix=prefix, sample=['A','B'])
    params:
        cells_meta = "cells_meta",
        #mt = pspace_mt.instance#(['mincells']) #pspace_cln.instance.get("n_clones")
    script: join(ROOT_DIR, "src/donor_filter_mgatk.py")
    # run:
    #     for p in params.mt:
    #         print(p)
    #         print(params.mt[p])
   # shell: "python src/tests/test_params.py {params.mt} && touch {output}"



# rule a:
#
#     output:
#          f"output/tests/pipeline/{prefix}/results/merged/{pspace_mt.wildcard_pattern}/dendrograms/af_dendro.ipynb"
#     params:
#         mt = pspace_mt.instance#(['mincells']) #pspace_cln.instance.get("n_clones")
#     # run:
#     #     for p in params.mt:
#     #         print(p)
#     #         print(params.mt[p])
#    # shell: "python src/tests/test_params.py {params.mt} && touch {output}"
#
