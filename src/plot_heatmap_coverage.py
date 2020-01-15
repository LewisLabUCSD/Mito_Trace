import pandas as pd
import glob
import os
import tqdm


def get_filter_CB():

    return

def create_heatmap_matrix():

    return


def fill_df_coverage(df, pileup_dir, is_par=False):
    not_done = []
    for ind in tqdm(df.index):
        f = glob.glob(os.path.join(pileup_dir,"CB_" + ind + ".coverage.txt"))
        if len(f) == 0:
            print("Not here")
            not_done.append(ind)
        else:
            curr = pd.read_csv(f[0], header=None)
            for _, val in curr.iterrows():
                df.loc[ind, val[0]] = val[2]
    print(f"Number of missing files: {len(not_done)}")
    return df