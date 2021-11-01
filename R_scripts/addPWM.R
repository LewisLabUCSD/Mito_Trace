# Title     : TODO
# Objective : TODO
# Created by: ras
# Created on: 10/26/21
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
plan("multiprocess", workers = 16)
#options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
options(future.globals.maxSize = 8000 * 1024^2)

get.pwm <- function(integrated, genome, out_f=""){
  # extract position frequency matrices for the motifs
  pwm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
  )

  fa.file <- Rsamtools::FaFile(genome, index=sprintf("%s.fai", genome)) #,

  DefaultAssay(integrated) <- "ATAC"
  chrom.assay <- integrated[["ATAC"]]
  # add motif information
  chrom.assay <- AddMotifs(chrom.assay, genome = fa.file, pfm = pwm)
  integrated = SetAssayData(integrated, slot="motifs", Motifs(chrom.assay))
  integrated = RunChromVAR(integrated, genome=fa.file)

  ## Add motif chromvarnames
  DefaultAssay(se) <- "chromvar"
  chrom.var.names <- GetAssayData(se)
  DefaultAssay(se) <- "ATAC"
  motifs <- Motifs(se)
  row.names(chrom.var.names) <- sapply(row.names(chrom.var.names), function(x) {Motifs(se)@motif.names[[x]]})

  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  integrated[['chromvarnames']] <- CreateAssayObject(data = chrom.var.names)
  DefaultAssay(integrated) <- "chromvarnames"


  if (out_f != ""){
    writeRDS(file=out_f, obj=integrated)
  }
  return(integrated)
}

