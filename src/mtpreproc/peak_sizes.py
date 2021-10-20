import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from src.utils.parse_config import read_config_file
from os.path import join


def get_peak_sizes(peaks_f, out_f):
    peaks = pd.read_csv(peaks_f, sep='\t', header=None)
    peaks.columns = ["Chr", "Start", "End", "Cell", "Count"]
    peaks["Size"] = peaks["End"]-peaks["Start"]
    sns.distplot(peaks["Size"])
    plt.savefig(out_f)
    return


# def wrap_peak_sizes(config_f):
# #     config = read_config_file(config_f)
# #     for i in config["indir"]:
# #         curr_cfg = read_config_file(config["indir"][i])
# #         curr_samples = pd.read_csv(curr_cfg["samples_meta"], sep='\t')
# #         for c in curr_samples.index:
# #             get_peak_sizes(curr_samples.loc[c,], out_f=)
# #     return
