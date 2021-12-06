## This is a module that gets called later and has multiple steps that depend on each other


rule pre:

rule a:
    input:
    output:
        "{outdir}/{prefix}/foo.txt",
        "{outdir}/{prefix}/foo.csv"
    shell: "touch {output}"


rule b:
    input:
         a="{outdir}/{prefix}/foo.txt",  #"{prefix}/foo.txt",prefix=["hi","hey"]),
         b="{outdir}/{prefix}/foo.csv"
    output:
        "{outdir}/{prefix}/{sample}/output.txt",
    params:
        cells_meta = "cells_meta",
        #mt = pspace_mt.instance#(['mincells']) #pspace_cln.instance.get("n_clones")
    shell: "touch {output[0]}"