lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4/" #'/usr/lib/R/library'

.libPaths(lib.loc) # Move to front for loading libeary
install.packages("dplyr", lib=lib.loc)
install.packages("ggrastr", lib=lib.loc)
install.packages("SummarizedExperiment", lib=lib.loc)
install.packages("ggrastr", lib=lib.loc)
install.packages("BuenColors", lib=lib.loc)
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager", lib=lib.loc)
BiocManager::install("GenomicRanges")


library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggrastr, lib=lib.loc)
library(BuenColors)

indir <- "/data2/isshamie/mito_lineage/data/processed/mtscATAC/jan21_2021/MTblacklist/reanalysis/outs"
indir_mgatk <- "/data2/mito_lineage/data/processed/mttrace/CHIP_april08_2021_Croker/MTblacklist/Flt3l/MT/cellr_True/Flt3l_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/"
name <- "Flt3l"
outdir <- "/data2/mito_lineage/Analysis/annotation/output/data"

#peaks_pan_heme <-diffloop::bedToGRanges("../data/29August2017_EJCsamples_allReads_500bp.bed")
#peaks_cd34 <-diffloop::bedToGRanges("../data/granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.peaks.bed")
peaks_cd34 <-diffloop::bedToGRanges("/data2/mito_lineage/Analysis/annotation/data/granja_cd34/GSE129785_scATAC-Hematopoiesis-CD34.peaks.bed")


# function to get counts
getCountsFromFrags <- function(frag_gz_file,
                               peaks_gr,
                               barcodes){
  
  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(paste0(frag_gz_file)) %>% 
    data.frame() %>% filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
  # Get a denominator, per cell
  denom <- table(GenomicRanges::mcols(frags_valid)$V4)
  barcodes_found <- names(denom)
  
  # Get the overlaps with peaks
  ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)
  
  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)
  
  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes_found
  
  # Make a polished colData
  colData <- data.frame(
    sample = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  # Make sure that the SE can be correctly constructed
  stopifnot(all(colData$sample == colnames(m)))
  
  # Make summarized Experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m),
    colData = colData
  )
  return(SE)
}

# Import experiment to make an RDS of the SE and a basic qc plot
importExperiment <- function(exp, peaks, frip_threshold){
  
  qcdf <- fread(paste0(indir,"singlecell.csv"), header = TRUE, sep = ",") %>% 
    data.frame() %>% filter(cell_id != "None")
  bc <- as.character(qcdf$barcode)
  SE <- getCountsFromFrags(paste0(indir, "filtered_peak_bc_matrix", "matrix.mtx"), peaks, bc) #fragments.tsv.gz"), peaks, bc)
  df <- data.frame(colData(SE))
  
  if(FALSE){
    df$density <- get_density(log10(df$depth + 1), df$FRIP*100)
    p1 <- ggplot(df %>% arrange(density), aes(x = log10(depth + 1),  y = FRIP*100, color = density)) +
      geom_point_rast(size = 0.1) + scale_color_gradientn(colors = jdb_palette("brewer_celsius")) +
      pretty_plot() + L_border() + labs(x = "log10 # Fragments", y = "FRIP") +
      theme(legend.position = "bottom") + ggtitle(exp) +
      geom_vline(xintercept = 3.5, linetype = 2) +
      geom_hline(yintercept = 40, linetype = 2) +
      scale_y_continuous(limits = c(0,100))
  }
  # Import mtDNA DF
  cov_mtDNA <- data.frame(colData(readRDS(paste0(indir_mgatk, exp, ".variant.rds"))))
  vec <- as.numeric(cov_mtDNA$depth); names(vec) <- as.character(cov_mtDNA$sample)
  colData(SE)$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  df$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  df$keep <- log10(df$depth) >= 3 & df$FRIP >= frip_threshold & df$mtDNAcoverage >= 20
  
  SE2 <- SE[, df$keep]
  saveRDS(SE2, file = paste0(outdir, exp, ".rds"))
  write.table(df, file = paste0(outdir, exp, ".barcode_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  print(dim(SE2))
}

# Use a much less stringent threshold for CD34 since we are borrowing the peak set
# from granja et al for projections
#importExperiment("cd34_jan21", peaks_cd34, 0.25)
importExperiment("Flt3l", peaks_cd34, 0.25)
#importExperiment("CD34_H8", peaks_cd34, 0.25)
#importExperiment("CD34_G10", peaks_cd34, 0.25)
#importExperiment("PBMC_H9", peaks_pan_heme, 0.6)
#importExperiment("PBMC_H10", peaks_pan_heme, 0.6)

