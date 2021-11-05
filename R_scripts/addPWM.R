# # Title     : TODO
# # Objective : TODO
# # Created by: ras
# # Created on: 10/26/21
# library(motifmatchr)
# library(JASPAR2020)
# library(TFBSTools)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(future)
#
#
# get.pwm <- function(se, genome, out_f=""){
#   # plan("multiprocess", workers = 16)
#   # #options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
#   # options(future.globals.maxSize = 8000 * 1024^2)
#   # extract position frequency matrices for the motifs
#   pwm <- getMatrixSet(
#     x = JASPAR2020,
#     opts = list(species = 9606, all_versions = FALSE)
#   )
#
#   fa.file <- Rsamtools::FaFile(genome, index=sprintf("%s.fai", genome)) #,
#
#   DefaultAssay(se) <- "ATAC"
#   chrom.assay <- se[["ATAC"]]
#   # add motif information
#   chrom.assay <- AddMotifs(chrom.assay, genome = fa.file, pfm = pwm)
#   se = SetAssayData(se, slot="motifs", Motifs(chrom.assay))
#   se = RunChromVAR(se, genome=fa.file)
#
#   ## Add motif chromvarnames
#   DefaultAssay(se) <- "chromvar"
#   chrom.var.names <- GetAssayData(se)
#   DefaultAssay(se) <- "ATAC"
#   motifs <- Motifs(se)
#   row.names(chrom.var.names) <- sapply(row.names(chrom.var.names), function(x) {Motifs(se)@motif.names[[x]]})
#
#   # add the gene activity matrix to the Seurat object as a new assay and normalize it
#   se[['chromvarnames']] <- CreateAssayObject(data = chrom.var.names)
#   DefaultAssay(se) <- "chromvarnames"
#
#
#   if (out_f != ""){
#     writeRDS(file=out_f, obj=se)
#   }
#   return(se)
# }
#
