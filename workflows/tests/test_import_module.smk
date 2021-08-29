## This is a module that gets called later and has multiple steps that depend on each other

rule a:
    output: "{prefix}/foo.txt"
    shell: "touch {output}"


rule b:
    input:
         "{prefix}/foo.txt"  #"{prefix}/foo.txt",prefix=["hi","hey"]),
    output:
        "{prefix}/results/{sample}/output.txt",
    params:
        cells_meta = "cells_meta",
        #mt = pspace_mt.instance#(['mincells']) #pspace_cln.instance.get("n_clones")
    shell: "touch {output[0]}"