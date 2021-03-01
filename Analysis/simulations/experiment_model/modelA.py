import numpy as np
import numpy.random as random
import pandas as pd
from itertools import product
from scipy.stats import fisher_exact
# 1. Poisson model


def initial_population(clone_percent, num_cells):
    #print(clone_percent)
    pop = np.sort(np.random.choice(np.arange(len(clone_percent)), num_cells, p=[clone_percent[0], clone_percent[1]]))
    return pop


def grow_cells(curr_pop, clone_growth_poisson):
    clone_ids = set(curr_pop)
    for clone in clone_ids:
        curr_size = (curr_pop==clone).sum()
        # Sample using the clone poisson, and the shape is based on the number of cells of that clone in the initial populatoin
        new_clones = clone * np.ones(shape=[random.poisson(lam=clone_growth_poisson[clone], size=curr_size).sum(), ])
        curr_pop = np.append(curr_pop, new_clones)
    return np.sort(curr_pop)


def sample_cells(cells, n_sample):
    #print('cells',len(cells), cells)

    #print('n_sample', n_sample)
    return np.random.choice(cells, n_sample, replace=False)


def determine_fold_over_init(pre, post):
    return


def shuffle_labels():
    return


# def run_fold_change(init, control, flt3):
#     clone_ids = set(control)
#     fold = {}
#     for i in clone_ids:
#         base = (init==i).sum()
#         fold[i] = (((flt3==i).sum())/base)
#                   ((control==i).sum())
#
#     return


def run_fold_no_init(control, flt3, chip_id):
    clone_ids = set(control)
    fold = {}
    enrich = {}
    total_chip =  (chip_id==flt3).sum() + (chip_id==control).sum()
    for i in clone_ids:
        a = (flt3==i).sum()
        b = (control==i).sum()
        N = len(flt3)+len(control)
        fold[i] = ((flt3==i).sum())/((control==i).sum())
    return fold


def test_single_run():
    init_cells = 5000
    clone_percent = {0: 0.2, 1: 0.8}
    clone_poisson = {0: 1, 1: 0.5}
    clone_poisson_flt3 = {0: 1.5, 1: 1}
    seq_cells = 2500
    n_simulations = 1000
    #single_run()
    return


def single_run(p):
    #print('p', p)
    pop = initial_population(p["clone_percent"], num_cells=p["init_cells"])
    cytokine_pop = grow_cells(pop, p["clone_poisson"])
    flt_pop = grow_cells(pop, p["clone_poisson_flt3"])
    cytokine_sample = sample_cells(cytokine_pop, n_sample=p["n_seq_cells"])
    flt_sample = sample_cells(flt_pop, n_sample=p["n_seq_cells"])
    fold = run_fold_no_init(cytokine_sample, flt_sample, chip_id=0)
    return pop, cytokine_sample, flt_sample, fold


from joblib import Parallel, delayed
def single_run_iter(n, in_params):
    print('iteration #', n)
    for ind, val in in_params.iterrows():
        curr_out_d = os.path.join(out_d, "params_" + str(ind))
        if not os.path.exists(curr_out_d):
            os.mkdir(curr_out_d)
        results = pd.DataFrame(
            columns=["CHIP Init", "Normal Init", "CHIP Flt3 seq",
                     "Normal Flt3 seq", "CHIP Cytokine seq",
                     "Normal Cytokine seq", "CHIP fold",
                     "Normal fold"])  # , "Fisher p"])

        pop, cytokine_sample, flt_sample, fold = single_run(val)
        results.loc[ind] = [(pop == 0).sum(), (pop == 1).sum(),
                            (flt_sample == 0).sum(),
                            (flt_sample == 1).sum(),
                            (cytokine_sample == 0).sum(),
                            (cytokine_sample == 1).sum(), fold[0], fold[1]]
        #print(type(n))
        #results.index.name = "ID"
        results.to_csv(os.path.join(curr_out_d, "sim" + str(n) + ".csv"))

    return



def wrap_params(out_d, n_sim=1000):
    n_seq_cells = [1500, 2000, 2500]
    l_init_cells = [2000,5000]
    l_clone_percent = [{0: 0.2, 1: 0.8},
                     {0: 0.02, 1: 0.98},
                     {0: 0.1, 1: 0.9}]
    l_clone_poisson = [
                     {0: 1, 1: 1},
                     {0: 0.5, 1: 0.5},
                     {0: 1, 1: 0.5}]
    l_clone_poisson_flt3 = [{0:4, 1:1}, {0: 4, 1: 0.5},
                          {0: 1.5, 1: 1},
                          {0: 1, 1: 0.5},
                          {0: 1, 1: 1}]

    in_params = pd.DataFrame(list(
        product(l_init_cells, l_clone_percent,
                l_clone_poisson, l_clone_poisson_flt3, n_seq_cells)),
                             columns=["init_cells", "clone_percent",
                                      "clone_poisson",
                                      "clone_poisson_flt3",
                                      "n_seq_cells"])

    in_params.to_csv(os.path.join(out_d, "in_params.csv"))

    Parallel(n_jobs=24)(
        delayed(single_run_iter)(n, in_params) for n in range(n_sim))

    return


from glob import glob
import seaborn as sns


def get_in_params():
    return ["init_cells", "clone_percent","clone_poisson","clone_poisson_flt3","n_seq_cells"]


def get_results_columns():
    return ["CHIP Init", "Normal Init", "CHIP Flt3 seq",
                     "Normal Flt3 seq", "CHIP Cytokine seq",
                     "Normal Cytokine seq", "CHIP fold",
                     "Normal fold"]


def figs(in_d, out_d):
    #1 Heatmap, where x is lambda of CHIP flt3, and y is CHIP flt3/CHIP ctrl
    # colors are different n_seq (multiplex) and columns are % CHIP

    # ["init_cells", "clone_percent","clone_poisson","clone_poisson_flt3","n_seq_cells"]
    # df_concat["CHIP Poisson cytokine"] = val["clone_poisson"][0]
    # df_concat["CHIP Poisson flt3"] = val["clone_poisson_flt3"][0]
    # df_concat["CHIP Real fold"] = df_concat["CHIP Poisson flt3"] / \
    #                               df_concat["CHIP Poisson cytokine"]
    #
    metrics = convert_to_long(in_d, out_d, max_n=100)
    metrics["measured CHIP fold"] = metrics["CHIP fold"]
    g = sns.catplot(col="CHIP Poisson flt3", y="measured CHIP fold",
                    x="CHIP Real fold", hue="n_seq_cells",
                    row="CHIP Percent", data=metrics, kind='violin')
    g.savefig(os.path.join(out_d, "chip_foldchange.png"))


    g = sns.catplot(col="CHIP Poisson flt3", y="CHIP Flt3 seq",
                    x="CHIP Real fold", hue="n_seq_cells",
                    row="CHIP Percent", data=metrics, kind='violin')
    g.savefig(os.path.join(out_d, "chip_flt3_cells.png"))
    # g = sns.FacetGrid(data=metrics, col='CHIP Percent',
    #                   row="n_seq_cells",
    #                   legend_out=True, margin_titles=True)
    #
    # g.map(sns.violinplot, "CHIP Poisson flt3", "CHIP fold")
    #

    #g.add_legend()  # bbox_to_anchor=(1.15, 1))


    #
    # g = sns.FacetGrid(data=metrics, col='CHIP Percent',
    #                   row="n_seq_cells", legend_out=True,
    #                   margin_titles=True)
    # g.map(sns.violinplot, "CHIP Poisson flt3", "CHIP fold")
    # g.add_legend()  # bbox_to_anchor=(1.15, 1))
    # g.savefig(os.path.join(out_d, "chip_fold_.png"))
    return


def convert_to_long(in_d, out_d, max_n=10):
    in_params = pd.read_csv(os.path.join(in_d, "in_params.csv"))
    # Create long dataframe, where extra columns are for the in_params (and their ID), and another column for the iter number
    #full_results = pd.DataFrame(columns=)
    all_dfs = []
    for ind, val in in_params.iterrows():
        curr_dfs = []

        files = glob(os.path.join(in_d, "params_" + str(ind), "*"))
        for curr_f in files[:min(len(files), max_n)]:
            print('curr_f', curr_f)
            curr_df = pd.read_csv(curr_f)
            curr_df["param ID"] = ind
            curr_dfs.append(curr_df)
        df_concat = pd.concat(curr_dfs, ignore_index=True)
        for name, v in  val.iteritems():
            df_concat[name] = v

        #["init_cells", "clone_percent","clone_poisson","clone_poisson_flt3","n_seq_cells"]
        print(val)
        print(val["clone_poisson"][0])
        df_concat["CHIP Poisson cytokine"] = float(eval(val["clone_poisson"])[0])
        df_concat["CHIP Poisson flt3"] = float(eval(val["clone_poisson_flt3"])[0])
        df_concat["CHIP Real fold"] = df_concat["CHIP Poisson flt3"]/df_concat["CHIP Poisson cytokine"]
        df_concat["CHIP Percent"] = eval(val["clone_percent"])[0]
        all_dfs.append(df_concat)
    all_dfs = pd.concat(all_dfs, ignore_index=True)
    all_dfs.to_csv(os.path.join(out_d, "results_long.csv"))
    return all_dfs



def expand_wrap(f_in, added_params_df):
    """ Add new params if the folder has been used already

    :param f_in:
    :param added_params_df:
    :return:
    """
    return



import os
out_d = "data/modelA/"
#wrap_params(out_d, n_sim=100)
figs(out_d, out_d)

