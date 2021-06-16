import snakemake
#workdir: "/data2/mito_lineage"
report: "../docs/workflow/template.rst"
#configfile: "parameters/tests/test_report.yaml"
from os.path import dirname, join
rule all:
    input:
        expand("{res}/01/p1_{p1}__p2_{p2}/02/p3_{p3}.txt", res=config["prefix"],
               p1=config["p1"], p2=config["p2"], p3=config['p3']),

rule t01:
    params:
        p1=lambda wc: wc.p1,
        p2=lambda wc: wc.p2
    output: "{res}/01/p1_{p1}__p2_{p2}/out1.txt"
    shell: "touch {output}"


rule t02:
    input:
        rules.t01.output
    params:
        p3=lambda wc: wc.p3,
        p4='p4'
    output:
        report("{res}/01/p1_{p1}__p2_{p2}/02/p3_{p3}.txt", category="Specific", caption="../docs/workflow/template2.rst")
    shell: "touch {output} "


rule t03_nothing:
    output: report(expand("{res}/03/nothing.txt", res=config['prefix']), category="General")
    shell: "touch {output} "
# # This is to be run after
# rule report:
#     input:
#         expand("{{res}}/01/p1_{p1}__p2_{p2}/02/p3_{p3}.txt",
#                p1=config["p1"], p2=config["p2"], p3=config['p3'])
#     output: "{res}/test_report.html"
#     run: snakemake.snakemake(snakefile="workflows/test_report.smk", report=output[0])
#     #shell: "snakemake --report "
#
