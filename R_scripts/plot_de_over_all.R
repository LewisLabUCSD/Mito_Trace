#source("./seurat_utils.R")
plot.over.all.top3 <- function(se, all.names.sig, curr.outdir, name="allClusters"){
  all.names.sig <- unique(all.names.sig)
  
  # Run de plots
  de.plots(se, all.names.sig, curr.outdir, curr.name=paste0(name, "combinedDonors"), max.size=15)
  
  dot <- DotPlot(se, scale=FALSE,
                 features = head(all.names.sig,20), 
                 cluster.idents=T) + RotatedAxis() + ggtitle(name)
  dot
  ggsave(file.path(curr.outdir, paste0(name, "combinedDonors.top3de.png")))
  ggsave(file.path(curr.outdir, paste0(name, "combinedDonors.top3de.pdf")))
}

plot.over.all.pvalsMean <- function(se, all.pvals, curr.outdir, name){
  ## b) averaged top p-values (only average ones that were DE)
  sig.all.ordered <- sort(rowMeans(-log10(all.pvals),na.rm=T), decreasing=T)
  dot <- DotPlot(se, scale=F,
                 features = names(head(sig.all.ordered,20)), 
                 cluster.idents=T) + RotatedAxis() + ggtitle(name)
  dot
  ggsave(file.path(curr.outdir, paste0(name,"combinedDonors.pvalsOrdered.scaleF.png")))
  ggsave(file.path(curr.outdir, paste0(name,"combinedDonors.pvalsOrdered.scaleF.pdf")))
}
plot.over.all.pvalsMeanIncNA <- function(se, all.pvals, curr.outdir, name){
  ## b) averaged top p-values (only average ones that were DE)
  all.pvals.nafill <- all.pvals
  all.pvals.nafill[,] <- -log10(zoo::na.fill(all.pvals,1))
  
  sig.all.ordered <- sort(rowMeans(all.pvals.nafill), decreasing=T)
  dot <- DotPlot(se, scale=F,
                 features = head(names(sig.all.ordered),20), 
                 cluster.idents=F) + RotatedAxis() + ggtitle(name)
  
  ggsave(file.path(curr.outdir,
                   paste0(name,"combinedDonors.ovalsOrderedNA.scaleF.png")))
  ggsave(file.path(curr.outdir,
                   paste0(name, "combinedDonors.ovalsOrderedNA.scaleF.pdf")))
  dot
}

wrap.plot.over.all <- function(se, curr.outdir, all.pvals, all.names.sig, name){
  plot.over.all.pvalsMeanIncNA(se, all.pvals, curr.outdir, name)
  plot.over.all.pvalsMean(se, all.pvals, curr.outdir, name)
  plot.over.all.top3(se, all.names.sig, curr.outdir, name)
}