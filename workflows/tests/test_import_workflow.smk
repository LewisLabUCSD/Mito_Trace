from snakemake.utils import min_version
min_version("6.0")

## This is a module that gets called later and has multiple steps that depend on each other
rule all:
    input:
        # expand("output/tests/pipeline/{prefix}/long/{sample}/output.txt",
        #        prefix=["hi","hey"], sample=['A','B']),
        expand("output/tests/pipeline/{prefix}/long_{long}/{sample}/output.D.txt",
               prefix=["hi","hey"], sample=['A','B'], long=['X','Y'])

# rule c:
#     input: "output/tests/pipeline/{prefix}/foo.txt"
#     output:
#         "output/tests/pipeline/{prefix}/long_{long}/{sample}/output.txt"
#     shell: "touch {output}"
rule c:
    output: "output/tests/pipeline/{prefix}/foo.txt"
    shell: "touch {output}"


module test_import_module:
    snakefile: "./test_import_module.smk"

use rule * from test_import_module
use rule b from test_import_module as b with:
    output: "output/tests/pipeline/{prefix}/long_{long}/{sample}/output.txt"


rule d:
    input: "output/tests/pipeline/{prefix}/long_{long}/{sample}/output.txt"
    output:
        "output/tests/pipeline/{prefix}/long_{long}/{sample}/output.D.txt"
    shell: "touch {output}"


# rule d:
#     input:
#          "output/tests/pipeline/{prefix}/long/foo.txt"  #"output/tests/pipeline/{prefix}/foo.txt",prefix=["hi","hey"]),
#     output:
#         "output/tests/pipeline/{prefix}/long/{sample}/output.D.txt"
#     params:
#         cells_meta = "cells_meta",
#         #mt = pspace_mt.instance#(['mincells']) #pspace_cln.instance.get("n_clones")
#     shell: "touch {output[0]}"