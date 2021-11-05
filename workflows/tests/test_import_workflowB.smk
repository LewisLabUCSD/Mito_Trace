from snakemake.utils import min_version
min_version("6.0")

## This is a module that gets called later and has multiple steps that depend on each other
rule all:
    input:
        expand("output/tests/pipelineB/long_{long}/{prefix}/{sample}/output.D.txt",
               prefix=["hi","hey"], sample=['A','B'], long=['X','Y'])

# rule a_long:
#     output:
#         "output/tests/pipelineB/long_{long}/{prefix}/{prefix}.txt",
#         "output/tests/pipelineB/long_{long}/{prefix}/{prefix}.csv"
#     shell: "touch {output}"


module test_import_moduleB:
    snakefile: "./test_import_moduleB.smk"

use rule * from test_import_moduleB


# use rule a from test_import_moduleB with:
#     output:
#         "output/tests/pipelineB/long_{long}/{prefix}/foo.txt",
#         "output/tests/pipelineB/long_{long}/{prefix}/foo.csv"


use rule b from test_import_moduleB with:
    input:
        #test_import_moduleB.rules
        a="output/tests/pipelineB/long_{long}/{prefix}/foo.txt",
        b="output/tests/pipelineB/long_{long}/{prefix}/foo.csv"
    output:
        "output/tests/pipelineB/long_{long}/{prefix}/{sample}/output.txt"


rule d:
    input: "output/tests/pipelineB/long_{long}/{prefix}/{sample}/output.txt"
    output:
        "output/tests/pipelineB/long_{long}/{prefix}/{sample}/output.D.txt"
    shell: "touch {output}"


# rule d:
#     input:
#          "output/tests/pipelineB/{prefix}/long/foo.txt"  #"output/tests/pipelineB/{prefix}/foo.txt",prefix=["hi","hey"]),
#     output:
#         "output/tests/pipelineB/{prefix}/long/{sample}/output.D.txt"
#     params:
#         cells_meta = "cells_meta",
#         #mt = pspace_mt.instance#(['mincells']) #pspace_cln.instance.get("n_clones")
#     shell: "touch {output[0]}"