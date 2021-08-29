# Title     : TODO
# Objective : TODO
# Created by: ras
# Created on: 2021-08-27
library(GenomicRanges)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(data.table)
library(magrittr)

library(future)
plan()
#plan("multiprocess", workers = workers)
options(future.globals.maxSize = 8000 * 1024^2)


read.peaks <- function(exp, cellr_in){
    print('here')
    print(file.path(cellr_in, exp, "outs", "filtered_peak_bc_matrix", "peaks.bed"))
    peaks <- read.table(
      file = file.path(cellr_in, exp, "outs", "filtered_peak_bc_matrix", "peaks.bed"),
      col.names = c("chr", "start", "end")
    )
    # convert to genomic ranges
    gr <- makeGRangesFromDataFrame(peaks)
    return(gr)
}


create_frag <- function(exp, cellr_in, combined.peaks, multiplex_cells=NULL, samples=NULL, sample_names=NULL, cores=2){

    if(is.null(multiplex_cells)){
        barcode_path <- file.path(cellr_in, exp, "outs", "filtered_peak_bc_matrix", "barcodes.tsv")
        barcodes <- readr::read_tsv(barcode_path, col_names = F) # %>% tidyr::unite(barcode)
        barcodes <- as.data.frame(barcodes) %>%  tibble::column_to_rownames(var="X1") %>% tibble::add_column(proj=exp)
    } else {
        for (i in 1:length(samples)){
            if (samples[[i]]==exp){
                cond <- sample_names[[i]]
        } }
        print("donor barcodes")
        print('cond')
        print(cond)

        barcodes <- multiplex_cells[multiplex_cells$condition == cond , c('raw.ID', "ID")] %>% tibble::add_column(proj=exp)
        rownames(barcodes) <- barcodes[,"raw.ID"]
        print(head(barcodes))
    }

    frag_file <- file.path(cellr_in, exp, "outs", "fragments.tsv.gz")
    # quantify multiome peaks in the scATAC-seq dataset
    frags.curr <- CreateFragmentObject(path = frag_file,cells= rownames(barcodes))

    ## Quantify peaks
    curr.counts <- FeatureMatrix(
      fragments = frags.curr,
      features = combined.peaks,
      cells = rownames(barcodes),
      process_n = cores
    )

    ## Create the objects
    curr_assay <- CreateChromatinAssay(curr.counts, fragments = frags.curr)
    curr <- CreateSeuratObject(curr_assay, assay = "ATAC", project=exp, meta.data=barcodes)
    curr <- BinarizeCounts(curr)
    return(curr)
}

integration <- function(combined){
    combined <- RunTFIDF(combined)
    combined <- RunSVD(combined)
    combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
    DimPlot(combined, group.by = "proj", pt.size = 0.1)
    # find integration anchors
    integration.anchors <- FindIntegrationAnchors(
      object.list = allSE,
      anchor.features = rownames(allSE[1]),
      reduction = "rlsi",
      dims = 2:30
    )

    # integrate LSI embeddings
    integrated <- IntegrateEmbeddings(
      anchorset = integration.anchors,
      reductions = combined[["lsi"]],
      new.reduction.name = "integrated_lsi",
      dims.to.integrate = 1:30
    )

    # create a new UMAP using the integrated embeddings
    integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

    return(integrated)
}
createMerged <- function(cellr_in, samples, sample_names, outdir, multiplex_cells=NULL, cores=2){
    gr.full <- c(sapply(samples, read.peaks, cellr_in=cellr_in, USE.NAMES=F))
    gr.full.c <- gr.full[[1]]
    if(length(gr.full)>1){
        for (i in 2:length(gr.full)){
            print("peak")
            print(i)
          gr.full.c <- c(gr.full.c, gr.full[[i]])
        }
    }
    combined.peaks <- reduce(x = c(gr.full.c))

    # Filter out bad peaks based on length
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

    allSE <- sapply(samples, create_frag, cellr_in=cellr_in,
                    combined.peaks=combined.peaks, multiplex_cells=multiplex_cells,
                    cores=cores, samples=samples, sample_names=sample_names)

    # Merge
    # merge all datasets, adding a cell ID to make sure cell names are unique
    combined <- merge(
      x = allSE[[1]],
      y = unlist(allSE[2:length(allSE)],use.names=FALSE), #allSE[2:length(allSE)],
      add.cell.ids = sample_names
    )
    combined[["ATAC"]]
    combined <- FindTopFeatures(combined, min.cutoff = 20)

    #combined <- subset(combined, nCount_peaks > 2000 & nCount_peaks < 30000)
    # Save merged
    saveRDS(combined, file.path(outdir, paste0("allSamples.merged.rds")))

    # save integrated
    integrated <- integrate(combined)
    saveRDS(integrated, file.path(outdir, paste0("allSamples.integrated.rds")))
    # Plot:
    p1 <- DimPlot(combined, group.by = "dataset")
    p2 <- DimPlot(integrated, group.by = "dataset")

    ggsave(((p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))), file = file.path(outdir, "allSamples.integrated.png"))
}



# Input info
cellr_in <- "/data2/isshamie/mito_lineage/data/processed/mtscATAC/jan21_2021/MTblacklist"
samples <- "J2,P2"
sample_names <- "Flt3l,Control"

# Saving
outdir <- "/data/isshamie/mito_lineage/output/annotation/cd34norm/MTblacklist/mergedSamples" #"/data2/mito_lineage/Analysis/annotation/output/data/"

# Parameters
#nTop = 25000
cores = 16

samples <- unlist(strsplit(samples, ",")[[1]])
sample_names <- unlist(strsplit(sample_names, ","))
print('samples')
print(samples)
print('outdir')
print(outdir)
multiplex_cells <- "/data2/mito_lineage/data/processed/mttrace/jan21_2021/MTblacklist/merged/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/vireoIn/multiplex/cells_meta.tsv"
multiplex_cells <- read.table(file=multiplex_cells, sep='\t', header=T)
#multiplex_cells <- lapply(samples, function(x) multiplex_cells[multiplex_cells$condition == x, 'raw.ID'])
print('multiplex cells')
print(head(multiplex_cells))
#multiplex_conditions <- sapply(samples, func(x) multiplex_cells["raw ID",condition=x] )

createMerged(cellr_in, samples, sample_names, outdir, multiplex_cells, cores)
