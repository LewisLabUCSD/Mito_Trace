The DE used is from the cell-by-TF normalized matrix, grouped by the largest clones (representing at least 30% of the top clones)
The following code was used in Seurat v4 to generate the results. The Benjamini Hochberg method was then used to adjust for the p-value::
    a = pairs[1,i]
    b = pairs[2,i]
    da <- FindMarkers(
      object = se.filt,
      ident.1 = a,
      ident.2 = b,
      only.pos = FALSE,
      mean.fxn = rowMeans,
      logfc.threshold = logfc_thresh,
      min.pct = min_pct,
      #latent.vars=latent.vars,
      fc.name = "avg_diff")

    da$p_val_adj_BH <- stats::p.adjust(da$p_val, method = "BH", n = length(da$p_val))
    write.csv(da,
              file=file.path(curr.outdir, paste0("clones_",a,"__", b,".DE.TF.csv")), quote=F)

From Seurat::
    p_val_adj: Adjusted p-value, based on bonferroni correction using all genes in the dataset
    min.pct: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1

    logfc.threshold
    Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.

P adjusted is bonferroni::
    >>> _, bonf,_,_ = multipletests(x["p_val"].values, method='bonferroni', is_sorted=False)
    >>> x["bonf"] = bonf
    >>> x.head()
        Unnamed: 0     p_val  avg_diff  pct.1  pct.2  p_val_adj  p_val_adj_BH      bonf
    0    JUN::JUNB  0.000025 -0.597928  0.330  0.474   0.015694       0.01254  0.015694
    1        FOSL2  0.000048 -0.654800  0.311  0.434   0.030182       0.01254  0.030182
    2  FOSL2::JUND  0.000093 -0.634327  0.324  0.425   0.058862       0.01254  0.058862
    3    FOS::JUND  0.000113 -0.634929  0.324  0.421   0.071385       0.01254  0.071385
    4   FOSL1::JUN  0.000124 -0.633695  0.333  0.439   0.078485       0.01254  0.078485
