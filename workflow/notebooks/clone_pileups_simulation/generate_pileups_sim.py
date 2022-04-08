import pandas as pd
import numpy as np
from collections import Counter

def generate_nonspec_vars(cell_ser, cell_nm, seq_err=0.001,
                          chars=["A", "C", "G", "T"]):
    #    print(cell_ser)
    reads = cell_ser["counts"]
    ref = cell_ser["ref"]
    ref_id = cell_ser["ref_id"]
    seq_err_counts = np.random.binomial(reads, seq_err)
    oth_nts = list(set(chars) - set([ref]))
    out = {}
    if seq_err_counts > 0:
        # Change the nt to the non alt allele
        out = Counter(np.random.choice(oth_nts, size=seq_err_counts,
                                       replace=True))
    else:
        out = {x: 0 for x in oth_nts}
    out[ref] = reads - seq_err_counts

    out["pos"] = cell_ser["pos"]
    out["cell"] = cell_nm

    return pd.Series(out)


def cell_nonspec_variants(curr_counts, cell_name):
    curr_cell_pile = curr_counts.apply(generate_nonspec_vars,
                                       args=(cell_name,),
                                       axis=1).reset_index().rename(
        {"index": "ref_id"}, axis=1)
    curr_cell_pile = curr_cell_pile.drop("ref_id", axis=1).melt(
        id_vars=["pos", "cell"], value_name="counts", var_name="nt")
    curr_cell_pile = (
    curr_cell_pile.loc[curr_cell_pile["counts"] != 0]).dropna()
    return curr_cell_pile


def generate_specific_vars(d_v_ser, cell_nm, don_var_lim, seq_err=0.001,
                           chars=["A", "C", "G", "T"]):
    # out = {}
    reads = d_v_ser["counts"]
    ref = d_v_ser["ref"]
    alt = d_v_ser["alt"]
    pos = d_v_ser["pos"]
    ref_id = d_v_ser["ref_id"]

    curr_af = np.random.uniform(don_var_lim[0], don_var_lim[
        1])  # generate using uniform distribution
    curr_af_counts = int(np.floor(curr_af * reads))
    # curr_af_counts = int(np.floor(curr_af*counts.loc[ref_id, "counts" ]))

    seq_err_counts = np.random.binomial(reads, seq_err)

    oth_nts = list(set(chars) - set([alt, ref]))
    if seq_err_counts > 0:
        # Change the nt to the non alt allele
        out = Counter(np.random.choice(oth_nts, size=seq_err_counts,
                                       replace=True))
    else:
        out = {x: 0 for x in oth_nts}

    # Add in the reference and alt counts!
    out[alt] = reads - seq_err_counts
    out[ref] = max(reads - curr_af_counts - seq_err_counts, 0)
    out["pos"] = d_v_ser["pos"]
    out["cell"] = cell_nm

    # Add the reference counts

    return pd.Series(out)


def cell_donor_variants(curr_cell_counts, cell_name, don_var_lim):
    """ Generates cell's donor variant counts based on the donor variant limits.
    """
    curr_cell_pile = curr_cell_counts.apply(generate_specific_vars,
                                            args=(
                                            cell_name, don_var_lim),
                                            axis=1).reset_index().rename(
        {"index": "ref_id"}, axis=1)
    #     print('curr_cell_pile')
    #     print(curr_cell_pile.head())
    # print(curr_cell_pile)
    curr_cell_pile = curr_cell_pile.drop("ref_id", axis=1).melt(
        id_vars=["pos", "cell"], value_name="counts", var_name="nt")
    curr_cell_pile = (
    curr_cell_pile.loc[curr_cell_pile["counts"] != 0]).dropna()
    return curr_cell_pile


def cell_clone_variants(curr_cell_counts, cell_name, clone_var_lim):
    """ Generates cell's counts based on the clone limits.
    """
    curr_cell_pile = curr_cell_counts.apply(generate_specific_vars,
                                            args=(
                                            cell_name, clone_var_lim),
                                            axis=1).reset_index().rename(
        {"index": "ref_id"}, axis=1)
    curr_cell_pile = curr_cell_pile.drop("ref_id", axis=1).melt(
        id_vars=["pos", "cell"], value_name="counts", var_name="nt")
    curr_cell_pile = (
    curr_cell_pile.loc[curr_cell_pile["counts"] != 0]).dropna()
    return curr_cell_pile


def cell_variants(cell_ser, ref_df, don_vars_df, clone_vars_df,
                  seq_err=0.001, depth_lim=(2, 10), strand_bin=0.5,
                  don_var_lim=(0.8, 1), clone_var_lim=(0.1, 0.5)):
    # cell_pileups = {}
    curr_clone = cell_ser["clone"]
    curr_don = cell_ser["donor"]
    cell_name = cell_ser.name

    # Get the donor and clone variants
    curr_don_vars = don_vars_df.loc[
        don_vars_df["donor"] == curr_don].set_index("ref_id")
    curr_cl_vars = clone_vars_df.loc[
        clone_vars_df["clone"] == curr_clone].set_index("ref_id")

    # Generate counts at each position
    counts = pd.DataFrame({"counts": np.floor(
        2 ** (np.random.randint(2, 10, size=ref_df.shape[0]))),
                           "ref": ref_df["ref"], "pos": ref_df["pos"],
                           "ref_id": ref_df.index, "cell": cell_name},
                          index=ref_df.index)
    counts["counts"] = counts["counts"].astype(int)

    # Non-donor and non-clone vars
    # print('non-spec counts')
    #     print('don vars', curr_don_vars.index)
    #     print('clone vars', curr_cl_vars.index)
    non_spec_counts = counts.drop(curr_don_vars.index).drop(
        curr_cl_vars.index)
    non_spec_pileup = cell_nonspec_variants(non_spec_counts, cell_name)

    # print('don spec')
    # First construct donor variants
    cell_donSpec_counts = counts.loc[curr_don_vars.index]
    cell_donSpec_counts["alt"] = curr_don_vars["alt"]
    cell_donSpec_pileup = cell_donor_variants(cell_donSpec_counts,
                                              cell_name, don_var_lim)

    # print('clone spec')
    # Construct clone variants
    cell_cloneSpec_counts = counts.loc[curr_cl_vars.index]
    cell_cloneSpec_counts["alt"] = curr_cl_vars["alt"]
    cell_cloneSpec_pileup = cell_clone_variants(cell_cloneSpec_counts,
                                                cell_name,
                                                clone_var_lim)

    #     for d_v in curr_don_vars:
    #         curr_don_pileup = {}
    #         curr_ref = curr_don_vars.loc[d_v, "ref"]
    #         curr_alt = curr_don_vars.loc[d_v, "alt"]
    #         curr_pos = curr_don_vars.loc[d_v, "pos"]
    #         curr_ref_id = curr_don_vars.loc[d_v, "ref_id"]
    #         curr_af = np.random.uniform(clone_var_lim[0],clone_var_lim[1])
    #         curr_af_counts = int(np.floor(curr_af*counts.loc[curr_ref_id, "counts" ]))
    #         seq_err_counts = np.random.binomial(counts.loc[curr_ref_id, "counts"],seq_err)
    #         curr_af_counts -= seq_err_counts
    #         if len(seq_error_counts)>0:
    #             # Change the nt to the non alt allele
    #             seq_err_counts = Counter(np.random.choice(set(chars)-set([curr_alt]), size=seq_err_counts,replace=True))

    #     print(cell_donSpec_counts.head())
    #     print(cell_cloneSpec_counts.head())
    #     print(non_spec.head())

    # concat donor-, clone-, and non-specific pileup dfs
    out_df = pd.concat(
        [cell_donSpec_pileup, cell_cloneSpec_pileup, non_spec_pileup],
        axis=0)
    out_df["donor"] = cell_ser["donor"]
    out_df["condition"] = cell_ser["condition"]
    return out_df

