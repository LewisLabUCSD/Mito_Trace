library(Signac)
library(Seurat)

filtCells <- function(se, min_peak_region_fragments=10,
                      max_peak_region_fragments=15000,
                     min_pct_reads_in_peaks=15,
                     max_blacklist_ratio=0.05,
                     max_nucleosome_signal=4,
                     min_TSS_enrichment=2){
    se <- subset(
      x = se,
      subset = peak_region_fragments > min_peak_region_fragments &
        peak_region_fragments < max_peak_region_fragments &
        pct_reads_in_peaks > min_pct_reads_in_peaks &
        blacklist_ratio < max_blacklist_ratio &
        nucleosome_signal < max_nucleosome_signal  &
        TSS.enrichment > min_TSS_enrichment
    )
        return(se)
}


embed.atac <- function(se, outdir, lsi_start_comp=2, reduction='lsi'){
# Binarize and run LSI
    se <- BinarizeCounts(se)
    se <- RunTFIDF(se)
    se <- RunSVD(se)
    se <- RunUMAP(se, dims = lsi_start_comp:50, reduction = reduction)
    DimPlot(se, group.by = "proj", pt.size = 0.1)
    pDepthCorr <- DepthCor(se, reduction=reduction)
    ggsave(file.path(outdir,"integrated.depthCor.png"), plot=pDepthCorr, dpi=300)
    pDepthCorr
        #integrated <- RunUMAP(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
    se <- FindNeighbors(object = se, reduction = 'lsi', dims = lsi_start_comp:30)
    se <- FindClusters(object = se, verbose = FALSE, algorithm = 3)
}


vPlot <- function(se){
      vPlot <- VlnPlot(
      object = se,
      features = c('pct_reads_in_peaks', 'peak_region_fragments',
                   'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
      pt.size = 0.1,
      ncol = 5
    )
    vPlot <- vPlot +    # Create grid of plots with title
             plot_annotation(title = se$orig.ident[[1]]) &
             theme(plot.title = element_text(hjust = 0.5, size=15))
    #print(vPlot)
    return(vPlot)
}