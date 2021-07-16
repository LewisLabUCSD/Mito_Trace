"""
File names for each rule, along with file loaders. Borrows from class templateDS which has Analysis-Pipeline-Protocol, where
Analysis is made of a list of Pipelines, made of Protocols. The Pipeline typically corresponds to a snakemake rule (or chain of rules),
and the Protocol is if there is subset of actions to take within the Pipeline (which could include mlflow).

Each pipeline overwrites the Loader class that takes as input the config, the pipeline, and vars to load.
"""
import pandas as pd



def paramspace():
    return


class Pipeline(object):
    def __init__(self, pipelines, config_f, samples=None, skeleteon_f=""):
        self.pipelines = pipelines
        self.config_f = config_f
        self.skeleton_f=skeleteon_f

    @private
    def check_set(self, var):
        if var == "config"
            return self.config == None
        if self.var is None:
            return False
        else:
            return True
    def set_samples(self):
        return

class mttrace(object):
    def __init__(self, config, type, samples=None):
        self.config=config

    def pipelines(self):
        return