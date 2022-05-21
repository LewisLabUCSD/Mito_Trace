# Title     : TODO
# Objective : TODO
# Created by: ras
# Created on: 5/18/22
#in_rds = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/annotation/gff_A2_black/mergedSamples/allSamples.integrated.rds"
se_f = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_3/gff_A2_black/annotation_clones/SE.rds"
labels_meta = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_3/gff_A2_black/annotation_clones/se_cells_meta_labels.tsv"
outdir = "/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/merged/MT/cellr_True/numread_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/mgatk/vireoIn/clones/variants_init/knn/kparam_3/gff_A2_black/annotation_clones/pseudotime"
order_f = ""
to_de = FALSE
#"/data/Mito_Trace/output/pipeline/v02/CHIP_b1/MTBlacklist_A2/data/annotation/gff_A2_black/mergedSamples/"
labels.meta <- read.table(labels_meta, sep="\t")



library(monocle3)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)

se <- readRDS(se_f)
se <- AddMetaData(se, labels.meta["cluster_labels"])
DefaultAssay(se) <- "ATAC"
se.cds <- as.cell_data_set(se)
se.cds <- cluster_cells(cds = se.cds, reduction_method = "UMAP")
cds_subset <- choose_cells(se.cds)

se.cds <- learn_graph(cds_subset, use_partition = TRUE)
# plot trajectories colored by pseudotime
plot_cells(
  cds = se.cds,
  show_trajectory_graph = TRUE, label_principal_points = TRUE
)
# order cells
# interactive more or not
if (order_f == ""){
   se.cds <- order_cells(se.cds, reduction_method = "UMAP", root_pr_nodes = c('Y_168','Y_157', 'Y_127'))
}else{
    se.cds <- order_cells(se.cds, reduction_method = "UMAP", root_cells = hsc)
}

# plot trajectories colored by pseudotime
plot_cells(
  cds = se.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
ggsave(file.path(outdir, "SE.pseudotime.trajectory.png"))

se <- AddMetaData(
  object = se,
  metadata = se.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)
FeaturePlot(se, c("pseudotime"), pt.size = 0.1) & scale_color_viridis_c()
ggsave(file.path(outdir, "SE.pseudotime.png"))
FeaturePlot(se, c("nCount_ATAC"), pt.size = 0.1) & scale_color_viridis_c()
ggsave(file.path(outdir, "SE.pseudotime.nCountPeaks.png"))


plot_cells(se.cds, color_cells_by="cluster_labels", label_cell_groups=TRUE)

plot_cells(se.cds, color_cells_by="name")
ggsave(file.path(outdir, "pseudo.clone.png"))

#cds_subset = cluster_cells(cds_subset, resolution=1e-2)
plot_cells(se.cds, color_cells_by="donor")
ggsave(file.path(outdir, "pseudo.donor.png"))


saveRDS(se, file.path(outdir, "SE.pseudotime.rds"))

saveRDS(se.cds, file.path(outdir, "SE.cds.rds"))