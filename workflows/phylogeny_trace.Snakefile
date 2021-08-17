

ref_fa = config["ref_fa"]
samples = pd.read_table(config["samples"], dtype=str,sep=',').set_index(["sample_name"], drop=False)
num_cells = config['multiplex']["pseudo_multiplex"]["num_cells"]
is_prop = config['multiplex']["pseudo_multiplex"]["is_proportional"]
res = config["results"]
ft = config["filters"]


rule all:
    input:
         expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/cass/treeclust/cells_BC.csv",
                    results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
                       het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),
         expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/cass/cells_phylo.gpickle",
                    results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
                       het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),


 expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/cass/treeclust/cells_BC.csv",
            results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),
 expand("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/cass/cells_phylo.gpickle",
            results=res, min_cells=ft['min_cells'],min_reads=ft['min_reads'],topN=ft["topN"],het_thresh=ft['het_thresh'],min_het_cells=ft['min_het_cells'],
               het_count_thresh=ft['het_count_thresh'], bq_thresh=ft['bq_thresh'], n_clones=config['multiplex']["n_clone_list"]),


def get_donor():
    return


rule run_cass:
    input:
        get_donor
    output:
        report("cass/output/data/cells_phylo_network.png"),
        report("cass/output/data/cells_phylo_hyperparams.png"),
        net = "cass/output/data/cells_phylo.gpickle",
        note = "cass/output/data/cells_phylo.ipynb",
    params:
        note="src/cass/01_reconstruct_from_af_df_HYBRID_MISSINGVALUES.ipynb"

    shell: "papermill -p input {input} -p output {output.net} {params.note} {output.note} && jupyter nbconvert --to pdf {output}"



rule run_treeclust:
    input: rules.run_cass.output.net
         #{raw_folder}/{bam_f}"
    output:
        cells = "cass/treeclust/output/data/cells_BC.csv",
        clust = report("cass/treeclust/output/data/cluster_sizes.png"),
        network_clust = report("cass/treeclust/output/data/cluster_net.png"),
        params_compare = report("cass/treeclust/output/data/suppl_cluster_params_elboA.png")
    params:
        note = "src/cass/02_cluster_from_dendrogram.ipynb"
    shell: "papermill -p input {input} -p output {output.cells} {params.note} {output.note} && jupyter nbconvert --to pdf {output}"


rule enrichment:
    input:
        #rules.multiplex.output[0],
        rules.clones.output[0],
    output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/enrichment/.status"
    params:
        clones_indir = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        n_clones = config['multiplex']["n_clone_list"],
        script = join("src", "vireo", "lineage_enrichment.py"),
        samples=",".join(config["multiplex"]['samples'])
    shell: "python {params.script} {params.clones_indir} {params.OUTDIR} {params.n_clones} {params.samples}"



rule clones_plotAF:
    input:
        #clone="data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/lineages/lineage.ipynb",
        mult = rules.multiplex.output[0]#"data/{prefix}/chrM/pseudo/minC{mt_minC}_minAF{mt_minAF}/numC{num_cells}_isprop{is_prop}/multiplex.ipynb"
    output:
        "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/dendrograms/dendrograms.ipynb",
        #report("{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/nfigures/donor{n}_dendrogram.png")
        report(expand("{{results}}/data/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/dendrograms/figures/donor{n}_dendrogram.png",
                n=np.arange(config["multiplex"]["N_DONORS"]),
                n_clones=config['multiplex']["n_clone_list"]))
    params:
        #INDIR = lambda wildcards, input: dirname(input[0]),
        INDIR = lambda wildcards, input: dirname(input.mult),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS=config["multiplex"]["N_DONORS"],
        sample_names= ','.join(samples["sample_name"].values), # make it as a list
        notebook=join("src", "vireo", "3_MT_Donors_Dendrogram.ipynb"),
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} {params.notebook} {output[0]}"


rule clones_type_variants:
    input: rules.clones.output[0]
    params:
        notebook=join("src", "vireo", join("6_MT_Clones_variantTypes.ipynb")),
        INDIR = lambda wildcards, input: dirname(input[0]),
        OUTDIR = lambda wildcards, output: dirname(output[0]),
        N_DONORS = config["multiplex"]["N_DONORS"],
        sample_names = ','.join(config['multiplex']["samples"]), # make it as a list
        n_clones = lambda wildcards: wildcards.n_clones, #config['multiplex']["n_clone_list"],#lambda wildcards: wildcards.n_clones,
        var_thresh=0.001,
        vars_to_plot=10
    #output: "{results}/merged/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/variants/variants.ipynb",
    output: "{results}/{sample}/clones/filters/minC{min_cells}_minR{min_reads}_topN{topN}_hetT{het_thresh}_hetC{min_het_cells}_hetCount{het_count_thresh}_bq{bq_thresh}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb"
    #output: report(lambda wildcards: expand("{{results}}/merged/filters/minC{{min_cells}}_minR{{min_reads}}_topN{{topN}}_hetT{{het_thresh}}_hetC{{min_het_cells}}_hetCount{{het_count_thresh}}_bq{{bq_thresh}}/filter_mgatk/vireoIn/clones/n_clones_{n_clones}/variants.ipynb", n_clones=config['multiplex']["n_clone_list"]))
    shell: "papermill -p INDIR {params.INDIR} -p OUTDIR {params.OUTDIR}  -p n_clones {params.n_clones} -p sample_names {params.sample_names} -p N_DONORS {params.N_DONORS} -p var_thresh {params.var_thresh} -p vars_to_plot {params.vars_to_plot} {params.notebook} {output} && jupyter nbconvert --to pdf {output}"

#rule: