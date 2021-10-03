#!/usr/bin/env Rscript

##!/usr/bin/Rscript --vanilla

# suppressMessages(suppressWarnings(library(tools)))
# suppressMessages(suppressWarnings(library(Matrix)))
# suppressMessages(suppressWarnings(library(SummarizedExperiment)))
# suppressMessages(suppressWarnings(library(GenomicRanges)))
# suppressMessages(suppressWarnings(library(data.table)))

suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(SummarizedExperiment)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(data.table)))
#library(BuenColors)
"%ni%" <- Negate("%in%")

call_mutations_mgatk <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10){
  
  # Determinie key coverage statistics every which way
  cov <- assays(SE)[["coverage"]]
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))
  
  # Process mutation for one alternate letter
  process_letter <- function(letter){
    print(letter)
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(SE))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    nucleotide <- paste0(ref_allele, ">", letter)[boo]
    position_filt <- pos[boo]
    
    # Single cell functions
    getMutMatrix <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_fw  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]]) / cov_fw)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_rev  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_rev")]]) / cov_rev)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    # Bulk functions
    getBulk <- function(letter){
      vec <- (Matrix::rowSums(assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / Matrix::rowSums(cov))[boo]
      return(vec)
    }
    rowVars <- function(x, ...) {
      Matrix::rowSums((x - Matrix::rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
    }
    
    update_missing_w_zero <- function(vec){
      ifelse(is.na(vec)  | is.nan(vec), 0, vec)
    }
    # Set up correlation per non-zero mutation based on the strands
    dt <- merge(data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_fw")]][boo,])), 
                data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_rev")]][boo,])), 
                by.x = c("i", "j"), by.y = c("i", "j"), 
                all = TRUE)[x.x >0 | x.y >0]
    dt$x.x <- update_missing_w_zero(dt$x.x)
    dt$x.y <- update_missing_w_zero(dt$x.y)
    
    dt2 <- data.table(variant = variant_name[dt[[1]]],
                      cell_idx = dt[[2]], 
                      forward = dt[[3]],
                      reverse = dt[[4]])
    rm(dt)
    cor_dt <- dt2[, .(cor = cor(c(forward), c(reverse), method = "pearson", use = "pairwise.complete")), by = list(variant)]
    
    # Put in vector for convenience
    cor_vec_val <- cor_dt$cor
    names(cor_vec_val) <- as.character(cor_dt$variant )
    
    # Compute the single-cell data
    mat <- getMutMatrix(letter)
    mmat <- sparseMatrix(
      i = c(summary(mat)$i,dim(mat)[1]),
      j = c(summary(mat)$j,dim(mat)[2]),
      x = c(update_missing_w_zero(summary(mat)$x), 0)
    )
    
    # Compute bulk statistics
    mean = update_missing_w_zero(getBulk(letter))
    
    # Stablize variances by replacing low coverage cells with mean
    if(stabilize_variance){
      
      # Get indices of cell/variants where the coverage is low and pull the mean for that variant
      idx_mat <- which(data.matrix(cov[boo,] < low_coverage_threshold), arr.ind = TRUE)
      idx_mat_mean <- mean[idx_mat[,1]]
      
      # Now, make sparse matrices for quick conversion
      ones <- 1 - sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = 1
      )
      
      means_mat <- sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = c(idx_mat_mean, 0)
      )
      
      mmat2 <- mmat * ones + means_mat
      variance = rowVars(mmat2)
      rm(mmat2); rm(ones); rm(means_mat); rm(idx_mat); rm(idx_mat_mean)
      
    } else {
      variance = rowVars(mmat)
    }
    
    detected <- (assays(SE)[[paste0(letter, "_counts_fw")]][boo,] >= 2) + (assays(SE)[[paste0(letter, "_counts_rev")]][boo,] >=2 )
    
    # Compute per-mutation summary statistics
    var_summary_df <- data.frame(
      position = position_filt,
      nucleotide = nucleotide, 
      variant = variant_name,
      vmr = variance/(mean + 0.00000000001),
      mean = round(mean,7),
      variance = round(variance,7),
      n_cells_conf_detected = Matrix::rowSums(detected == 2),
      n_cells_over_5 = Matrix::rowSums(mmat >= 0.05), 
      n_cells_over_10 = Matrix::rowSums(mmat >= 0.10),
      n_cells_over_20 = Matrix::rowSums(mmat >= 0.20),
      strand_correlation = cor_vec_val[variant_name],
      mean_coverage = Matrix::rowMeans(cov)[boo], 
      stringsAsFactors = FALSE, row.names = variant_name
    )
    se_new <- SummarizedExperiment(
      rowData = var_summary_df, 
      colData = colData(SE), 
      assays = list(allele_frequency = mmat, coverage = cov[boo,])
    )
    return(se_new)
  }
  
  return(SummarizedExperiment::rbind(process_letter("A"), process_letter("C"), process_letter("G"), process_letter("T")))
  
}




plot_mutations_qc <- function(mut_se, f_save, is_df = FALSE) {
  
  if (is_df == FALSE) { 
    misc_df <- data.frame(rowData(mut_se))  
  } else {
    misc_df <- mut_se
  }
  
  # filter_df <- misc_df %>%  filter(n_cells_conf_detected >= n_cells_thresh & strand_correlation >strand_correlation_thresh & log10(vmr) > log_vmr_thresh)
  # dim(filter_df)
  # filter_df # Verify that 8202 and 8344 are there
  #
  # Make the standard variant calling plot
  p1 <- ggplot(misc_df, aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
    geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
    labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
    #L_border() + pretty_plot(fontsize = 8) + ## Not in R 3.5
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0.65, linetype = 2) +
    geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
  cowplot::ggsave2(p1, file = paste(gsub('.png','', f_save),".png", sep=""), width = 1.7, height = 1.7)
}



# i/o
####################
#-----------------
# hard-coded i/o
#-----------------
#SE_f <- "/data2/mito_lineage/data/processed/mttrace/jan21_2021/P2/MT/cellr_True/P2_200/filters/minC100_minR50_topN0_hetT0.01_hetC10_hetCount5_bq30/filter_mgatk/P2.rds"#"/data2/mito_lineage/data/processed/mttrace/jan21_2021/P2/MT/cellr_True/P2_200/filter_mgatk/P2.rds"
# SE_f <- "/data2/mito_lineage/data/processed/mttrace/2020_11_18/PBMC_J/mapq_0/cellr_True/PBMC_J_200/PBMC_J.rds"
# SE_f <- "/data2/mito_lineage/data/processed/mttrace/jan21_2021/P2/MT/cellr_True/P2_200/mgatk/P2.rds" #"P2.rds" 
#out_f <- "tmp_variant_old.rds"
#low_coverage_threshold <- 0

#-----------------
# Command line i/o
#-----------------
# args <- commandArgs(trailingOnly = TRUE)
# #print(sessionInfo())
# if (length(args) == 3) {
#   SE_f <- args[1]
#   low_coverage_threshold <- args[2]
#   n_cells_thresh <- args[3] #2
#   print(low_coverage_threshold)
#   ####################
# 
# 
#   
#   SE <- readRDS(SE_f)
#   strand_correlation_thresh <- 0.65
#   log_vmr_thresh <- -2
#   
#   
#   out_SE <- gsub('.rds', paste("_lowC", low_coverage_threshold, "_cellT", n_cells_thresh, sep = ""), SE_f)
#   print('out_SE')
#   print(out_SE)
#   
#   if (!(grepl('.rds', SE_f))) {
#     print('Not an .rds file! Not running')
#   } else {
#     # Call variants
#     
#     mut_se <- call_mutations_mgatk(SE, low_coverage_threshold=low_coverage_threshold)
#     print('mut_se')
#     print(head(mut_se))
#     
#     saveRDS(mut_se, file = gsub('.rds', '.variant.rds', SE_f))
#   
#     
#     misc_df <- data.frame(rowData(mut_se))
#     filter_df <- misc_df %>%  filter(n_cells_conf_detected >= n_cells_thresh & strand_correlation > strand_correlation_thresh & log10(vmr) > log_vmr_thresh)
#     #write.table(as.data.frame(as.matrix(assay(mut_se, 2))), file = paste(out_SE, ".coverage.tsv", sep = ""), sep='\t')
#     cov <- as.data.frame(as.matrix(assay(mut_se, 2)))
#     cov -> cov %>% filter(row.names(cov) %in% row.names(filter_df))
#     write.table(cov, file = paste(out_SE, ".coverage.tsv", sep = ""), sep='\t')
#     af <- as.data.frame(as.matrix(assay(mut_se, 1)))
#     af  -> af %>% filter(row.names(af) %in% row.names(filter_df))
#     write.table(af, file = paste(out_SE, ".af.tsv", sep = ""), sep='\t')
#     write.table(filter_df, file = paste(out_SE, ".af.mgatk.tsv", sep = ""), sep='\t')
#     plot_mutations_qc(mut_se , f_save = paste(out_SE,'.variantQC.png', sep = ""))
#   }
# } else {
#   print("Args not correct: SE_f; low_coverage_threshold; n_cells_thresh")
# }
  

# # Find TF1 cells
# rbind(read.table("../output/data1_meta.tsv", header = TRUE), 
#       read.table("../output/data6_meta.tsv", header = TRUE)) %>% filter(assign == "GM11906") %>%
#   filter(mean_cov > 50) -> GM_cells_df
# 
# SE <- SE[,as.character(GM_cells_df$cell_id)]


