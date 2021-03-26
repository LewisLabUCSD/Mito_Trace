# Clone-Cell type relationship Writeup
- Once we get the clones, and we correlate them with cell types, we want to see if there's any relationship between the clones and the cell types
- This will be similar to Lareau et al 2020, and Lin et al 2021.

Current Report:
---

## A. Data:
- CD34+ cells 
- Condition A: cytokine cocktail; 3 days
- Condition B: Flt3l + cytokine cocktail; 3 days


## B. Clustering cells on scATAC-seq peaks to enable cell type calling
[Control](./clusters_P2.png)

[Flt3l](./clusters_J2.png)

![Merged](./clusters_reanalysis.png)

![Merged Conditions](./reanalysis_aggr_conditions.png)


### Marker DC TFs from Lin

![Id2](./id2.png)

![Irf8](./reanalysis_irf8.png)


### TFs enriched in clusters
This csv file provides statistics for enriched clusters:
[TF](./Graph_Based_Genes_allsignificant.csv)


### Plotting markers as heatmap across cells 
[Immune markers](./immune_markers_peak_results_merged_conditions.pdf)

Note the clusters for the cells are from the 10x clusters, which
was over the entire genome, not just near these markers)

These were the markers used. Some might be incorrect (sca1). 
- cKit:KIT,
- Sca1: LY6E, #Ly6a,LY6K LY6E LY6H ??
- CD11c: ITGAX,
- CD150: SLAMF1,
- CD34: CD34, 
- CD16/32:FCGR3A,
- CD45.1: PTPRC, 
- CD45.2:PTPRC,
- CD48: SLAMF2, # Other SLAMFs 
- IL7Ra:IL7R
- CD11b:ITGAM


## C. Multiplexing Donors
[Donors in 'cell type' embedding](./cells_merged_lin_and_peak_donors.tsne.subplots.png)

## D.  Calling clones on the cells

See Google slides for now.

## E. Grouping cells by both cell-types and Clones

![Clone Lineage 20 clones each donor](./cells_merged_lin_and_peak_nclones20.overlap_percent_normClone.png)

[Clone Lineage 100 clones each donor](./cells_merged_lin_and_peak_nclones100.overlap_percent_normClone.png)

Note that one donor did not have any confident clones called (p>0.9 a cell is assigned to a clone) when n-clones is 100.

## F. Inspired by Lin et al Fig2. , dim-red cell-type:clone cell count matrix highlights lineages that are certain found certain 'fate-cluster', so we can find some lineages

- Donors are separated
- There are still dynamics can be seen between lineages and cell types

['Fate-Clusters'](./clone_fate_scanpy_separateConditions.pdf)

- UMAP of Clone-Lineage 20 clones each donor.
    - Each lineage (flt3l and wt are separated) are their own datapoint, with number of cell-type clusters being the number of features. All donors included

['Fate-Clusters' separating out the donors](./clone_fate_scanpy_separateConditions_SplitDonors.pdf)

## To Do:

- A: The 10x seems to not mix flt3l and WT very well.
    - Can try a few different parameters from 10x or a UMAP method , but the umap shouldn't matter too much
    - Another alternative is to use cell similarity score based on prior scATAC-seq values.
    - Can also use marker genes for this
- C: Still need to nail down the clones
- E. This looks very interesting and promising, but need good cell type annotations and flt3l-wt mix
