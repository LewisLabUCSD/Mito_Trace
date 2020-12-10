
plot_mutations_qc <- function(mut_se, f_save, strand_correlation_thresh = 0.65, n_cells_thresh = 5, log_vmr_thresh=-2) {
  misc_df <- data.frame(rowData(mut_se))
  filter_df <- misc_df %>%  filter(n_cells_conf_detected >= n_cells_thresh & strand_correlation >strand_correlation_thresh & log10(vmr) > log_vmr_thresh)
  dim(filter_df)
  filter_df # Verify that 8202 and 8344 are there
  
  # Make the standard variant calling plot
  p1 <- ggplot(filter_df, aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
    geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
    labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
    #L_border() + pretty_plot(fontsize = 8) + ## Not in R 3.5
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0.65, linetype = 2) +
    geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
  cowplot::ggsave2(p1, file = paste(f_save,".png", sep=""), width = 1.7, height = 1.7)
}
