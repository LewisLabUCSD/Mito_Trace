"""File names using params"""
from os.path import join
import pandas as pd
import itertools

def create_f(p, rule, sample=None):
    #print('p', p)
    print('rule', rule)
    # if rule == "main":
    #     return f"{p.outdir}/{p.prefix}/data/"
    if rule == "cellrbc_sample":
        return f"{sample}/MT/cellr_{p.cellrbc}/numread_{p.numreadsfilter}/"
    elif rule == "cellrbc":
        return f"merged/MT/cellr_{p.cellrbc}/numread_{p.numreadsfilter}/"
    elif rule == "create_filters":
        return f"filters/minC{p.mincells}_minR{p.minreads}_topN{p.topN}_hetT{p.hetthresh}_hetC{p.minhetcells}_hetCount{p.hetcountthresh}_bq{p.bqthresh}/"
    elif rule == "mgatk":
        return "mgatk/vireoIn"
    elif rule == "multiplex":
        return "multiplex"

    elif rule == "clones":
        if p.method == "knn":
            return f"clones/variants_{p.variants}/{p.method}/kparam_{p.resolution}"
        elif p.method == "vireo":
            return f"clones/variants_{p.variants}/{p.method}/nclones{p.nclonelist}"

    elif rule == "enrichment":
        return "enrichment/volcano_Fisher_foldNorm.png"

    elif rule == "annotation_clones":
        return f"annotation_clones/DE_large/minPct_{p.min_pct}__logThresh_{p.logfc_threshold}/cdf_thresh__{p.cdf_thresh}"
        #return f"{p['output']}/data/merged/MT/cellr_{p}/numread_{num_read}/filters/minC{mincells}_minR{minreads}_topN{topN}_hetT{hetthresh}_hetC{minhetcells}_hetCount{hetcountthresh}_bq{bqthresh}/mgatk/vireoIn/clones/variants_mgatkdonor/knn/kparam_{kparam}/enrichment/volcano_Fisher_foldNorm.png",
    else:
        raise ValueError(f"No rule {rule}")
    return


def wrap_create_f(params_ser, cfg, rule, sample=None):
    if rule == "main":
        return join(cfg["outdir"], "pipeline", cfg["prefix"], "data")
    elif rule == "cellrbc":
       # print('main run', wrap_create_f(params_ser, cfg, "main"))
        #print('create cellrbc', create_f(params_ser, "cellrbc"))
        return join(wrap_create_f(params_ser, cfg, "main"),
                    create_f(params_ser, "cellrbc"))
    elif rule == "cellrbc_sample":
        return join(wrap_create_f(params_ser, cfg, "main"),
                    create_f(params_ser, "cellrbc_sample", sample=sample))
    elif rule == "create_filters":

        return join(wrap_create_f(params_ser, cfg, "cellrbc"),
                    create_f(params_ser, "create_filters"))

    elif rule == "mgatk":
        return join(wrap_create_f(params_ser,cfg, "create_filters"),
                    create_f(params_ser, "mgatk"))

    elif rule == "multiplex":
        return join(wrap_create_f(params_ser, cfg, "mgatk"),
                    create_f(params_ser, "multiplex"))

    elif rule == "clones":
        return join(wrap_create_f(params_ser, cfg, "mgatk"),
                    create_f(params_ser, "clones"))

    elif rule == "enrichment":
        return join(wrap_create_f(params_ser, cfg, "clones"),
                    create_f(params_ser, "enrichment"))

    elif rule == "annotation_clones":
        print('here')
        return join(wrap_create_f(params_ser, cfg, "clones"),
                    create_f(params_ser, "annotation_clones"))


def create_single_files(cfg, rule, sample=None):
    all_p = {}
    all_p['filt'] = cfg["filters"]["params"]
    all_p['mtpreproc'] = cfg["mtpreproc"]["params"]
    #print(cfg["clones"]["method"])
    all_p['mgatk'] = cfg["mgatk"]["params"]
    all_p['clones'] = {}
    all_p['clones']["method"] = cfg["clones"]["method"]
    all_p['clones']["variants"] = cfg["clones"]["variants"]
    for m in cfg["clones"]["method"]:
        for curr_p in cfg["clones"][m]["params"]:
            all_p['clones'][curr_p] = cfg["clones"][m]["params"][curr_p]
    all_p['annotation_clones'] = cfg["annotation_clones"]["params"]
    #methods_df = pd.DataFrame

    cols = []
    col_vals = []
    col_d = {}
    for p in all_p:
        col_d[p] = []
        for k in list(all_p[p].keys()):
            if k in cols:
                #print('new col', f"{k}_{p}")
                cols += [f"k_{p}"]
                col_d[p].append(f"{k}_{p}")
            else:
                cols += [k]
                col_d[p].append(k)
            col_vals.append(tuple(all_p[p][k]))
    print('cols', cols)
    print('col_vals', col_vals)
    print(len(cols))
    print(len(col_vals))
    print('iter')
    #print(list(itertools.product(*col_vals)))
    params_df = pd.DataFrame(list(itertools.product(*(col_vals))), columns=cols)
                            #columns=["mtpreproc","filt", "mgatk", "clones", "annotation_clones"])
    print("cfg", cfg)
    params_df["file"] = params_df.apply(wrap_create_f, cfg=cfg, rule=rule, sample=sample, axis=1)
    return params_df


# """
# File names for each rule, along with file loaders. Borrows from class templateDS which has Analysis-Pipeline-Protocol, where
# Analysis is made of a list of Pipelines, made of Protocols. The Pipeline typically corresponds to a snakemake rule (or chain of rules),
# and the Protocol is if there is subset of actions to take within the Pipeline (which could include mlflow).
#
# Each pipeline overwrites the Loader class that takes as input the config, the pipeline, and vars to load.
# """

# import pandas as pd
# from src.utils.data_io import import_from_dir
# from src.calculate_AF_by_cell import load_and_filter_nt_pileups
#
#
# class Pipeline(object):
#     def __init__(self, pipelines, config_f, samples=None, skeleteon_f=""):
#         self.pipelines = pipelines
#         self.config_f = config_f
#         self.skeleton_f=skeleteon_f
#
#     def check_set(self, var):
#         if var == "config":
#             return self.config == None
#         if self.var is None:
#             return False
#         else:
#             return True
#
#     def set_samples(self):
#         return
#
#
#
#
# class mttrace(object):
#     def __init__(self, num_read, samples=None, version="v1"):
#         if version == "v1":
#             self.vars = {"scPileup": {"suffix": f"{num_read}_",
#                                       "names":samples},
#                          }
#         else:
#             raise ValueError(f"Version {version} not available. Please check again")
#         return
#
#     def get_files(self, indir=".", vars=None, names=None):
#         curr_files = {}
#         if vars is None:
#             vars = self.vars
#         for v in vars:
#             curr_files[v] = import_from_dir(indir=indir, prefix=v.get("prefix", None), suffix=v.get("suffix", None),
#                                            f_list=names)
#         return
#
#     def load_files(self, indir, var_name, file_name):
#         if var_name == "scPileup":
#             for i in ["A", "C", "T", "G", "coverage"]:
#
#                 return load_and_filter_nt_pileups(indir, cell_inds, nt_pos,
#                                out_d=None, name=None, incl_cov=False,
#                                avg_cov=False)(file_name)
#         return
#
#
#
#
#
# PIPELINES = {
#     "mttrace": mttrace,
#     "multiplex": mttrace,
#     "cluster": mttrace,
#     "annotation": annotation,
#     "cluster__graphclust": mttrace
# }
#
#
