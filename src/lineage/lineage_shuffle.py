""" Script to calculate statistics for the fold-change of the enrichment.

"""


def get_fold(shuff, a, b):
    a_df = shuff[shuff["condition"]==a]
    b_df = shuff[shuff["condition"]==b]

    fold_df = {}
    return


def shuffle_lineage_clones(n, cells_meta, conditions):
    lineages = cells_meta["lineage"]

    n_lineage = max(cells_meta["lineage"].values)

    # Shuffle lineage assignment based on random index
    shuff = cells_meta.sample(frac=1).index

    shuff_clones = cells_meta.copy()
    shuff_clones["shuffle_lineage"] = cells_meta.loc[shuff, "lineage"]

    fold_df = get_fold(shuff_clones, conditions[0], conditions[1])
    return

