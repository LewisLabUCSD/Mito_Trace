#source("toRDS.R")
source("variant_calling.R")

low_coverage_threshold <- 10
n_cells_thresh <- 5
strand_correlation_thresh <- 0.65
log_vmr_thresh <- -2



"data/processed/mttrace/CHIP_april08_2021_Croker/MTblacklist/Flt3l/MT/cellr_True/Flt3l_200/filters/minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20/filter_mgatk/Flt3l.rds" 

################################################################################################
data_dir <- "/data2/mito_lineage/data/processed/mttrace/CHIP_april08_2021_Croker/MTblacklist"
samples <- c("Control", "Flt3l", "Input") 

indirs <- lapply(samples, function(x) file.path(data_dir, x, "MT", "cellr_True", paste0(x,"_200"), "filters", "minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20"))
names <- lapply(samples, function(x) file.path('mgatk',paste0(x)))
seurat_rds <- lapply(samples, function(x) file.path(data_dir, x, "MT", "cellr_True", paste0(x,"_200"), "filters", "minC10_minR50_topN0_hetT0.001_hetC10_hetCount5_bq20",'filter_mgatk', paste0(x,'.rds')))
# make directory if not exist 
lapply(seurat_rds, function(x) ifelse(!dir.exists(dirname(x)), dir.create(dirname(x)), FALSE))
all <- data.frame(samples=samples, indirs=unlist(indirs), names=unlist(names), rds=unlist(seurat_rds))
################################################################################################

wrap_to_seurat <- function(folder, name, is_strand){
  if (is_strand == "TRUE"){
    is_strand <- TRUE
  } else if (is_strand =="FALSE"){
    is_strand<-FALSE
  }
  ####################
  print(sessionInfo())
  
  SElist <- importMito(folder, is_strand)
  saveRDS(SElist[[1]], file = paste0(folder, "/", name, ".rds"))
  saveRDS(SElist[[2]], file = paste0(folder, "/", name, ".signac.rds"))
}


wrap_var <- function(SE_f) {
  SE <- readRDS(SE_f)
  if (!(grepl('.rds', SE_f))) {
    print('Not an .rds file! Not running')
  } else {
    # Call variants
    mut_se <- call_mutations_mgatk(SE, low_coverage_threshold=low_coverage_threshold)
    #misc_df <- data.frame(rowData(mut_se))
    #filter_df <- misc_df %>%  filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2)
    #dim(filter_df)
    print(mut_se)
    saveRDS(mut_se, file = gsub('.rds', '.variant.rds', SE_f))
    
    misc_df <- data.frame(rowData(mut_se))
    filter_df <- misc_df %>%  filter(n_cells_conf_detected >= n_cells_thresh & strand_correlation > strand_correlation_thresh & log10(vmr) > log_vmr_thresh)
    write.table(as.data.frame(as.matrix(assay(mut_se, 2))), file = gsub('.rds', ".coverage.tsv", SE_f), sep='\t')
    write.table(as.data.frame(as.matrix(assay(mut_se, 1))), file = gsub('.rds', "af.tsv", SE_f), sep='\t')
    write.table(filter_df, file = gsub('.rds', "af.mgatk.tsv", SE_f), sep='\t')
    plot_mutations_qc(mut_se , f_save = gsub('.rds', 'variantQC.png', SE_f))
  }
}

print("all")
print(all)

#apply(all, MARGIN=1, FUN=function(x) wrap_to_seurat(folder=x["indirs"],name=x["names"], is_strand = TRUE))
apply(all, MARGIN=1, FUN=function(x) wrap_var(x["rds"]))
