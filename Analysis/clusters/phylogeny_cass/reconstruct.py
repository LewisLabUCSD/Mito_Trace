import pandas as pd
import numpy as np
import cassiopeia as cas
import seaborn as sns
from os.path import join, exists
from os import makedirs, getcwd
from pandarallel import pandarallel
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import write_dot, graphviz_layout

pandarallel.initialize(nb_workers=32)


af_f = "/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/post/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/af_by_cell.tsv"
dp_f = "/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/post/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/af_by_cell.DP.tsv"
ad_f = "/data2/mito_lineage/data/processed/mttrace/TcellDupi_may17_2021/MTblacklist/post/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/af_by_cell.AD.tsv"
prefix= "TcellDupi_may17_2021/MTblacklist"
name = "post_cass_test"
outdir = "./output/data/"
af_thresh = 0.01
dp_thresh = 2

outdir = join(outdir, prefix)
if not exists(outdir):
    print(f"Making outdir {outdir} in folder {getcwd()}")
    makedirs(outdir)
    
af = pd.read_csv(af_f, index_col=0, sep='\t')
dp = pd.read_csv(dp_f, index_col=0, sep='\t').astype(int)
ad = pd.read_csv(ad_f, index_col=0, sep='\t').astype(int)

#Binarize
af = af.applymap(lambda x: 0 if x<af_thresh else 1)
af
#af[af>=af_thresh] = 1
#af = af.loc[dp<dp_thresh] = -1

def dp_where(x, dp, dp_thresh):
#    print('name', x.name)
    curr = dp.loc[x.name]<dp_thresh
#     print('curr')
#     print(curr)
    x.loc[curr] = -1 #"-"
    return x

#af = af.parallel_apply(dp_where, axis=1, args=(dp,))
# tmp0 = af.head(10).apply(dp_where, axis=1, args=(dp,0))

# tmp10 = af.head(10).apply(dp_where, axis=1, args=(dp,10))
# (tmp0==-1).sum(axis=1)
character_matrix = af.copy()
#character_matrix = af.parallel_apply(dp_where, axis=1, args=(dp,dp_thresh))
character_matrix


var_map = {val:f"r{i}" for i, val in enumerate(character_matrix.columns) }
character_matrix = character_matrix.rename(var_map, axis=1)

sample_mat = character_matrix.sample(n=100)
sample_mat = sample_mat.loc[:, ~((sample_mat==-1).all(axis=0))]
sample_mat = sample_mat.loc[~((sample_mat==-1).all(axis=1)), :]

var_map = {val:f"r{i}" for i, val in enumerate(sample_mat.columns) }
sample_mat = sample_mat.rename(var_map, axis=1)

sample_mat.to_csv(join(outdir, name+".sample100.tsv"),sep='\t')

sample_mat.shape
priors=None
if priors is not None:
    cas_tree = cas.data.CassiopeiaTree(character_matrix=sample_mat, priors=priors)
else:
    cas_tree = cas.data.CassiopeiaTree(character_matrix=sample_mat)
cas_tree.character_matrix.head(5)
# REINSTANTIATE the bottom and top solvers
vanilla_greedy = cas.solver.VanillaGreedySolver()
if priors is not None:
    ilp_solver = cas.solver.ILPSolver(convergence_time_limit=500, maximum_potential_graph_layer_size=500, weighted=True, seed=1234)
else:
    ilp_solver = cas.solver.ILPSolver(convergence_time_limit=500, maximum_potential_graph_layer_size=500, weighted=False, seed=1234)
    
# hybrid_solver = cas.solver.HybridSolver(top_solver=vanilla_greedy, bottom_solver=ilp_solver, cell_cutoff=40, threads=10)
# hybrid_solver.solve(cas_tree, logfile='hybrid.log')
hybrid_solver = cas.solver.HybridSolver(top_solver=vanilla_greedy, bottom_solver=ilp_solver, cell_cutoff=40, threads=10)
hybrid_solver.solve(cas_tree, logfile='example_hybrid.log')
