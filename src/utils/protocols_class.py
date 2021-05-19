
from src.utils.parse_config import read_config_file, check_required
from src.config import ROOT_DIR
from mplh.fig_utils import helper_save as hs

import mlflow

import os
from os.path import join
import pandas as pd
import logging



def set_artifact():
    ###TODO: add this into the class
    if not os.path.exists("output"):
        os.mkdir("output")
    for f in ["images", "plots", "data", "results"]:
        if not os.path.exists(join("output",f)):
            os.mkdir(join("output",f))



def helper_add_prefix_to_recursive_dict(d, in_dir):
    """ Recursive function to add a prefix to every value in the dict

    :param d:
    :param prefix:
    :return:
    """
    out = {}
    for key in d:
        if type(d[key]) is dict:
            out[key] = helper_add_prefix_to_recursive_dict(d[key], in_dir)
        elif type(d[key]) is not str:
            raise ValueError("Incorrect files key")
        else:
            out[key] = join(in_dir, d[key])
    print(out)
    return out


def test_helper_add_prefix_to_recursive_dict():
    assert(helper_add_prefix_to_recursive_dict({"h":"data"},
                                               in_dir="test/")=={"h":"test/data"})
    d2 = {"h":"data", "h2": {"h2.1":"hi.csv", "h2.2": "hey.csv"}}
    assert(helper_add_prefix_to_recursive_dict(d2, "/data2") == {
        "h": "/data2/data", "h2": {"h2.1": "/data2/hi.csv", "h2.2": "/data2/hey.csv"}})
    return


def helper_dict_flatten(d:dict, out_d:dict=None):
    """

    :param d:
    :param out_d:
    :return:
    """
    if out_d is None:
        out_d = {}
    for k in d:
        if type(d[k]) == dict:
            out_d = helper_dict_flatten(d[k], out_d=out_d)
        else:
            if k in out_d:
                print('overkap k', k)
                raise ValueError("There are overlapping keys")
            out_d[k] = d[k]

    return out_d

def test_helper_dict_flatten():
    assert(helper_dict_flatten({"h":"data"})=={"h":"data"})
    d2 = {"h":"data", "h2": {"h2.1":"hi.csv", "h2.2": "hey.csv"}}
    assert(helper_dict_flatten(d2) == {"h": "data",
                                       "h2.1": "hi.csv",
                                       "h2.2": "hey.csv"})
    return

test_helper_dict_flatten()
#test_helper_add_prefix_to_recursive_dict()

########
## Config
########
class Config:
    def __init__(self, params_d=None, file=None):
        print(f"Config file: {file}")

        if params_d is not None:
            self.params = params_d
            self.params_f = ""
        else:
            self.params_f = file
            self.params = read_config_file(file)

        # Verify certain necessary parameters
        check_required(self.params, ["Prefix", "Project_Directory",
                                     "Working_Directory"])
        self.prefix = self.params["Prefix"]
        # Setup default parameters
        # Setup file paths
        # Save some metadata and full paths
        if self.params["Project_Directory"] in [".", None, "Current", ""]:
            self.project_dir = ROOT_DIR
        else:
            self.project_dir = self.params["Project_Directory"]

        if self.params["Working_Directory"] in [".", None, "Current", ""]:
            self.curr_dir = os.getcwd()
        elif self.params["Project_Directory"] in ["ROOT_DIR", "Project"]:
            self.curr_dir = self.project_dir

        #self.absfiles = self.files
        return


    def create_full_path(self):
        return


class AnalysisConfig(Config):
    def __init__(self, file):
        super().__init__(file)
        # Verify certain necessary parameters
        check_required(self.params, ["Prefix"])

        # Setup default parameters
        self.protocols = self.params.get('protocols',[])

        # Setup file paths
        # Save some metadata and full paths
        self.project_dir = self.params["Project_Directory"]
        self.prefix = self.params["Prefix"]
        #self.absfiles = self.files


# class _ProtocolsConfig(Config):
#     def __init__(self, params=None, file=None):
#         super().__init__(params, file)
#         self.input_incfg = self.params["input"]
#         self.output_incfg = self.params["output"]
#         self.parameters = self.params["parameters"]
#         self.protocols_name = self.params["Protocol_Name"]
#
#         self.input_d = {}
#         self.output_d = {}
#
#     def initialize(self):
#         self.setup_input_files()
#         self.setup_output_files()
#         return
#
#     def setup_input_files(self):
#         curr = self.input_incfg
#         protocols = curr.get("protocol", [])
#         in_files = self.input_from_protocols(protocols)
#         return
#
#
#     def input_from_files(self):
#         return
#
#     def input_from_parameters(self):
#         return




########
## Protocol
########
class Protocol(object):
    def __init__(self, config_f=None, config=None,
                 pipeline_id=-1,
                 input_d=None,
                 output_d=None,
                 mlflow_tracking_uri="",
                 experiment_name=None):
        if config_f is not None:
            self.config = read_config_file(config_f)
        else:
            self.config = config
        self.config_f = config_f
        self.workdir = os.getcwd()
        self.vars = {}
        self.saved_f = {}
        if input_d is None:
            self.input_dict = {}
        else:
            self.input_dict = input_d
        if output_d is None:
            self.output_dict = {}
        else:
            self.output_dict = output_d
        self.pipeline_id = pipeline_id
        self.analysis = Analysis
        self.config_processed = False
        self.is_completed = False

        ## mlflow attributes
        self.mlflow_tracking_uri = mlflow_tracking_uri
        self.run_id = None
        self.experiment_name=experiment_name

        self.input_incfg = self.config["input"]
        self.output_incfg = self.config["output"]
        self.parameters = self.config["parameters"]

        self.samples_metadata_f = self.config.get("samples_meta", None)
        self.initialize_config()
        if self.samples_metadata_f is not None:
            samples_meta_is_project = self.parameters.get("samples_meta_is_project", True)
            if samples_meta_is_project:
                self.samples_metadata_f = join(self.project_dir, self.samples_metadata_f)
            self.metadata = pd.read_csv(self.samples_metadata_f)
            self.samples_meta_is_project = samples_meta_is_project
        return


    def run(self):
        self.initialize()
        self.load_input_files()
        self.init_mlflow()
        return


    def init_mlflow(self, mlflow_tracking_uri=None, run_id=None):
        """
        Sets the mlflow tracking uri, and adds run_id if provided.
        :param mlflow_tracking_uri:
        :param run_id:
        :return:
        """
        if mlflow_tracking_uri is not None:
            self.mlflow_tracking_uri = mlflow_tracking_uri

            if run_id is not None:
                self.run_id = run_id

        mlflow.set_tracking_uri(self.mlflow_tracking_uri)
        if self.experiment_name is not None:
            mlflow.set_experiment(self.experiment_name)
        return

    ############
    ## Initialize the config parameters and the input + output dicts
    ############
    def initialize(self):
        self.initialize_config()
        self.initialize_dicts()
        self.config_processed=True
        return

    def initialize_config(self):
        # if self.config is None:
        #     raise ValueError("Need to set config first in init")
        # self.input_incfg = self.config["input"]
        # self.output_incfg = self.config["output"]
        # self.parameters = self.config["parameters"]
        # print('config')
        # print(self.config)
        if self.config["Project_Directory"] in ["", None, ROOT_DIR]:
            self.project_dir = ROOT_DIR
        else:
            self.project_dir = self.config["Project_Directory"]

        if self.config["Project_Directory"] in [".", "", None]:
            self.working_directory = ""
        else:
            self.working_directory = self.config["Working_Directory"]

        self.protocols_name = self.config["Protocol_Name"]
        self.prefix = self.config["Prefix"]
        self.protocol_name = self.config["Protocol_Name"]
        self.workflow_file = self.config.get("workflow_file", "")
        return


    def initialize_dicts(self):
        #self.input_dict = self.config.setup_input_files()
        self.setup_input_files()
        self.setup_output_files()
        return


    #############
    ## Setup of the full input_dict and full output_dict,
    ## which are dicts with keys as file keys and values the actual
    ## file
    ############

    ## Input setup
    def setup_input_files(self):
        curr = self.input_incfg
        protocols = curr.get("protocols", {})
        in_files = self.input_from_protocols(protocols)
        for f in in_files:
            self.input_dict[f] = in_files[f]
        files = curr.get("files", {})
        in_files = self.input_from_files(files)
        for f in in_files:
            self.input_dict[f] = in_files[f]

        #self.input_dict = helper_dict_flatten(self.input_dict)
        return

    ## Output setup
    def setup_output_files(self):
        outdir = self.output_incfg.get("outdir", ".")
        is_proj = self.output_incfg.get("is_project", False)
        if is_proj:
            outdir = join(self.project_dir, outdir)
        self.output_dict = self.output_files(outdir=outdir,
                                             params=self.parameters)
        self.set_outdir(outdir, self.output_dict)
        return

    @staticmethod
    def set_outdir(outdir, output_dict=None):
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if output_dict is None:
            return
        else:
            for d in output_dict:
                if not os.path.exists(d):
                    os.mkdir(d)

    #############
    ## 2 ways to setup input: from explicit files and from protocols
    ############
    ## From files
    def input_from_files(self, files):
        """

        :param files:
        :return:
        """
        full_files = {}

        for f in files:
            is_proj = files[f].get("is_project", True)
            dir = files[f].get("dir", "")
            curr = ""
            if is_proj:
                curr = ROOT_DIR
            if dir:
                curr = join(curr, ROOT_DIR)
            full_files[f] = curr
        return full_files

    ## From protocols
    def input_from_protocols(self, protocols):
        """ Setup input files from the prior protocols

        :param protocols:
        :return:
        """
        pipeline_files = {}
        for p in protocols:
            curr = self.setup_protocol_files(protocols[p], p)
            print('protocol', p)
            print('curr', curr)

            for c in curr:
                if c in pipeline_files:
                    print(c)
                    raise ValueError("key already there")
                pipeline_files[c] = curr[c]
        return pipeline_files

    def setup_protocol_files(self, protocol, name):
        logging.log(logging.DEBUG, f"protoocol: {name}")
        is_project, dir, files, parameters = \
            [protocol.get('is_project', True),protocol.get('dir', ""),
             protocol.get('files', {}),
             protocol.get('parameters', {})]

        if is_project:
            dir = join(self.project_dir, dir)

        files = self.analysis.protocol_file_setup(protocol=name,
                                             outdir=dir,
                                             params=parameters,
                                             files=files
                                            )
        print("files", files)
        files = helper_add_prefix_to_recursive_dict(files, in_dir=dir)
        return files

    def load_input_files(self, var_keys=None, params=None,
                         reload=False, return_vars=True):
        """ Loads the vars keys from the input_dict attr using var_loader method.

        :param vars:
        :return:
        """
        if params is None:
            params=self.parameters

        var_keys = self.setup_var_keys(var_keys)
        input_keys = self.input_keys()
        for v in var_keys:
            if (not input_keys[v] in self.vars) or reload:
                self.vars[input_keys[v]] = self.var_loader(
                    input_keys[v], self.input_dict[v], params=params)
            elif reload:
                print(f"{v} already loaded")
        if return_vars:
            return self.vars
        return

    #########################################################
    ## Unique for each Protocol
    ########################################################
    @staticmethod
    def input_keys():
        return {}

    ###############
    @staticmethod
    def var_loader(var, fname, params=None):
        """

        :param var:
        :param fname:
        :param params:
        :return:

        if var == "af":
            return read_af(fname)
        elif var == "mgatk_af":
            return read_af(fname)
        else:
            return f"No loader for variable {var}"
        """
        return f"No loader for variable {var}"
        # cov = read_af(in_cov_d)



    ###############
    @staticmethod
    def output_files(outdir="", output_dict=None, params=None, name=""):
        """ This is unique for each protocol. If output_d is provided,
        it will overwrite those keys that are there.

        :param outdir: The outdirectory for each output ifle.
        :param output_dict: Dictionary with keys results, plots, data that
        point to different file names for certain output.
        :return:
        """
        curr_output_dict = {"results": {}, "data": {"test.txt"},
                                "plots": {}}
        # add directory information
        for d in ["results", "data", "plots"]:
            for f in curr_output_dict[d]:
                curr_output_dict[d][f] = join(outdir, curr_output_dict[d][f])

        # Prior dictionary assumes the paths are full
        if output_dict is None:
            return curr_output_dict
        for c in output_dict:
            if c not in ["results", "plots", "data"]:
                raise ValueError("Output keys are either results, plots, data")
            for v in output_dict[c]:
                curr_output_dict[c][v] = output_dict[c][v]
        return curr_output_dict

    ## Till Here
    ########################################################
    ########################################################
    def load_var(self, v):
        """
        Make specific cases for the variables. Additionally, tests
        if the variable is a proper type.
        :param v:
        :return:
        """
        if v == "df list":
            pd.read_csv(v)
        return

    def save(self, var, outtype,vars=None):
        join(outtype, self.vars[var].file)
        return

    def add_var(self, x, out_f, name):
        self.vars[name] = Var(out_f, x)
        return


    def setup_var_keys(self, var_keys):
        """ Gets the variable names and their associated files

        :param var_keys:
        :return:
        """
        if self.input_dict is None:
            raise AttributeError("input_dict is not set yet")
        if var_keys is None:
            var_keys = list(self.input_dict.keys())
        return var_keys

    def change_to_wkdir(self, workdir=None):
        if workdir is None:
            self.workdir = os.getcwd()
        else:
            if os.path.exists(workdir):
                self.workdir = workdir
                os.chdir(workdir)
            else:
                raise ValueError(f"{workdir} does not exist.")

    def setup_absolute_files(self):
        return




########
## Logging
########
class Log:
    def __init__(self, log_f=None):
        logging.basicConfig(filename=log_f, filemode='w', level='DEBUG',
                            format='%(name)s - %(levelname)s - %(message)s')
        logging.info(f"Using log file {log_f}")
        return


########
## protocol vars
########
class Var:
    def __init__(self, file=None, var=None):
        """ A variable that has a load and save, a value, and file/files associated with it.
        The load may be flexible based on different types of input.
        :param var:
        :param file:
        """
        print('file', file)
        self.file = file

        self.var = var
        self.saved_in_this_run=False
        self.loaded = False
        if file is None:
            self.is_file = False
        else:
            self.is_file=True
        return

    def save(self, f_save, **kwargs):
        return

    def load(self, file=None):
        """ Loads the variable, if the file is present

        :param file:
        :return:
        """
        pass
        if file is not None:
            self.file = file
        if self.file is None:
            raise ValueError("Please set file first")
        return


class Output(Var):
    def __init__(self, file, var):
        super().__init__(file, var)


class Report(Output):
    def __init__(self, file, var):
        super().__init__(file, var)
        return


class Plot(Output):
    def save(self, f_save, **kwargs):
        hs(f_save, **kwargs)
        return


class Data(Output):
    def __init__(self, file, var):
        super().__init__(file, var)
        self.outdir = "data/processed"
        return



########
## Analysis
########
class Analysis:
    def __init__(self, params_f, name, description=None):
        self.set_metadata_csv()
        self.description = description
        self.params_f = params_f
        self.config = AnalysisConfig(params_f)

        self.protocols = self.config.protocols
        self.pipeline = []
        self.params = ""
        self.project_dir = ""
        self.analysis_dir = ""
        self.name = name


        #self.params = self.load_params(params_f)
        return

    def run(self):
        for p in self.protocols:
            self.protocol_runner(p)
        return

    @staticmethod
    def protocol_file_setup(outdir="", protocol=None, params=None,
                                    files=None):
        """ Setup the output files for the protocols

        Setup a specific protocol if it is
        not None.
        This is done for each specific anlayses.
        An example is shown below.

        :param outdir: directory to use for the output
        files, which will be input into each protocol
        :param protocol: The name of the specific protoocl
        to use.
        :param params:
        :return:

        if params is None:
            params = {}
        loaders = {"MTExp" : MTExp.output_files(outdir=outdir,
                                                params=params)}
        if protocol is None:
            return loaders
        return loaders[protocol]
        """
        if protocol is None:
            return None
        return {}

    @staticmethod
    def protocol_runner(protocol=None):
        """ Setup the protocols to be used

        Run a specific protocol if it is
        not None.
        This is done for each specific anlayses.
        An example is shown below.

        :param protocol:
        :return:

        runners = {"MTExp": MTExp().run()}
        if protocol is None:
            return runners
        return runners[protocol]
        """
        if protocol is None:
            return None
        return {}
    def log_output(self):
        return

    def load_params(self, params_f):
        return

    ######
    ## Parameters setup methods
    ######
    def set_metadata_csv(self):
        self.metadata = pd.read_csv(self.config.params["samples_meta"])
        return


    def output_log(self):
        return

    ######
    ## helper methods
    ######
    def add_output(self, output):
        return

    def check_if_files_are_written(self):
        return

    def pipeline_add(self, params):
        self.pipeline.append(params)
        return

