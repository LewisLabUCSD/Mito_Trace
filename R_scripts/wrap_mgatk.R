#!/usr/bin/Rscript --vanilla
source("./R_scripts/toRDS.R")
source("./R_scripts/variant_calling.R")

low_coverage_threshold <- 10
n_cells_thresh <- 5
strand_correlation_thresh <- 0.65
log_vmr_thresh <- -2

wrap_to_seurat <- function(folder, name, is_strand){
  if (is_strand == "TRUE"){
    is_strand <- TRUE
  } else if (is_strand =="FALSE"){
    is_strand<-FALSE
  }
  ####################
  print(sessionInfo())
  SElist <- importMito(folder, is_strand)
  print(paste0(folder, "/", name, ".rds"))
  saveRDS(SElist[[1]], file = paste0(folder, "/", name, ".rds"))
  saveRDS(SElist[[2]], file = paste0(folder, "/", name, ".signac.rds"))
}


wrap_var <- function(SE_f) {
  if (!(grepl('.rds', SE_f))) {
    print('Not an .rds file! Not running')
  } else {
    SE <- readRDS(SE_f)
    # Call variants
    mut_se <- call_mutations_mgatk(SE, low_coverage_threshold=low_coverage_threshold)
    #misc_df <- data.frame(rowData(mut_se))
    #filter_df <- misc_df %>%  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2)
    #dim(filter_df)
    print('mut_se')
    print(mut_se)
    saveRDS(mut_se, file = gsub('.rds', '.variant.rds', SE_f))

    out_SE <- gsub('.rds','', SE_f) # paste("_lowC", low_coverage_threshold, "_cellT", n_cells_thresh, sep = ""), SE_f)
    misc_df <- data.frame(rowData(mut_se))
    filter_df <- misc_df %>%  filter(n_cells_conf_detected >= n_cells_thresh & strand_correlation > strand_correlation_thresh & log10(vmr) > log_vmr_thresh)
    #write.table(as.data.frame(as.matrix(assay(mut_se, 2))), file = paste(out_SE, ".coverage.tsv", sep = ""), sep='\t')
    cov <- as.data.frame(as.matrix(assay(mut_se, 2)))
    cov <- cov %>% filter(row.names(cov) %in% row.names(filter_df))
    write.table(cov, file = paste(out_SE, ".coverage.tsv", sep = ""), sep='\t')
    af <- as.data.frame(as.matrix(assay(mut_se, 1)))
    af  <- af %>% filter(row.names(af) %in% row.names(filter_df))
    write.table(af, file = paste(out_SE, ".af.tsv", sep = ""), sep='\t')
    write.table(filter_df, file = paste(out_SE, ".af.mgatk.tsv", sep = ""), sep='\t')
    plot_mutations_qc(mut_se , f_save = gsub('.rds', '.variantQC.png', SE_f))
  }
}

args <- commandArgs(trailingOnly = TRUE)
print(sessionInfo())
folder <- args[1]
name <- args[2]
is_strand <- args[3]
if (is_strand == "TRUE"){
  is_strand <- TRUE
} else{
  is_strand<-FALSE
}
wrap_to_seurat(folder, name, is_strand)
SE_f <- file.path(folder, paste0(name, ".rds"))
print("SE_f")
print(SE_f)
wrap_var(SE_f)
