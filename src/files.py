"""
File names for each rule, along with file loaders. Borrows from class templateDS which has Analysis-Pipeline-Protocol, where
Analysis is made of a list of Pipelines, made of Protocols. The Pipeline typically corresponds to a snakemake rule (or chain of rules),
and the Protocol is if there is subset of actions to take within the Pipeline (which could include mlflow).

Each pipeline overwrites the Loader class that takes as input the config, the pipeline, and vars to load.
"""
import pandas as pd
from src.utils.data_io import import_from_dir
from src.calculate_AF_by_cell import load_and_filter_nt_pileups


class Pipeline(object):
    def __init__(self, pipelines, config_f, samples=None, skeleteon_f=""):
        self.pipelines = pipelines
        self.config_f = config_f
        self.skeleton_f=skeleteon_f

    def check_set(self, var):
        if var == "config":
            return self.config == None
        if self.var is None:
            return False
        else:
            return True

    def set_samples(self):
        return




class mttrace(object):
    def __init__(self, num_read, samples=None, version="v1"):
        if version == "v1":
            self.vars = {"scPileup": {"suffix": f"{num_read}_",
                                      "names":samples},
                         }
        else:
            raise ValueError(f"Version {version} not available. Please check again")
        return

    def get_files(self, indir=".", vars=None, names=None):
        curr_files = {}
        if vars is None:
            vars = self.vars
        for v in vars:
            curr_files[v] = import_from_dir(indir=indir, prefix=v.get("prefix", None), suffix=v.get("suffix", None),
                                           f_list=names)
        return

    def load_files(self, indir, var_name, file_name):
        if var_name == "scPileup":
            for i in ["A", "C", "T", "G", "coverage"]:

                return load_and_filter_nt_pileups(indir, cell_inds, nt_pos,
                               out_d=None, name=None, incl_cov=False,
                               avg_cov=False)(file_name)
        return





PIPELINES = {
    "mttrace": mttrace,
    "multiplex": mttrace,
    "cluster": mttrace,
    "annotation": annotation,
    "cluster__graphclust": mttrace
}


