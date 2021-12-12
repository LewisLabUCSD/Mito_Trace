library(Signac)
library(Seurat)

filtCells <- function(se, min_peak_region_fragments=10,
                      max_peak_region_fragments=15000,
                     min_pct_reads_in_peaks=15,
                     max_blacklist_ratio=0.05,
                     max_nucleosome_signal=4,
                     min_TSS_enrichment=2){
    se <- subset(
      x = se,
      subset = peak_region_fragments > min_peak_region_fragments &
        peak_region_fragments < max_peak_region_fragments &
        pct_reads_in_peaks > min_pct_reads_in_peaks &
        blacklist_ratio < max_blacklist_ratio &
        nucleosome_signal < max_nucleosome_signal  &
        TSS.enrichment > min_TSS_enrichment
    )
        return(se)
}


embed.atac <- function(se, outdir, lsi_start_comp=2, reduction='lsi'){
# Binarize and run LSI
    se <- BinarizeCounts(se)
    se <- RunTFIDF(se)
    se <- RunSVD(se)
    se <- RunUMAP(se, dims = lsi_start_comp:50, reduction = reduction)
    DimPlot(se, group.by = "proj", pt.size = 0.1)
    pDepthCorr <- DepthCor(se, reduction=reduction)
    #ggsave(file.path(outdir,"integrated.depthCor.png"), plot=pDepthCorr, dpi=300)
    #pDepthCorr
        #integrated <- RunUMAP(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
    se <- FindNeighbors(object = se, reduction = 'lsi', dims = lsi_start_comp:30)
    se <- FindClusters(object = se, verbose = FALSE, algorithm = 3)
    return(se, pDepthCorr)
}


vPlot <- function(se){
      vPlot <- VlnPlot(
      object = se,
      features = c('pct_reads_in_peaks', 'peak_region_fragments',
                   'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
      pt.size = 0.1,
      ncol = 5
    )
    vPlot <- vPlot +    # Create grid of plots with title
             plot_annotation(title = se$orig.ident[[1]]) &
             theme(plot.title = element_text(hjust = 0.5, size=15))
    #print(vPlot)
    return(vPlot)
}


check.if.one.clone <- function(large.clones, init.large.clones, n_top_clones) {
    donor.n.clones <- large.clones %>% group_by(donor) %>% summarise(size=n())
    #' If one of the donors has only 1 donor, use the top n instead
    print(dim(donor.n.clones)[1])
    for (i in 1:dim(donor.n.clones)[1]) {
        curr.d <- donor.n.clones[i, "donor"]
        print(paste('i', i, "donor", curr.d))
        print(donor.n.clones[i, "size"])
        if (donor.n.clones[i, "size"] < n_top_clones) {
            print(paste('only 1 clone. Using top', n_top_clones, 'clones'))
            new.large <- init.large.clones %>% filter(donor==as.character(curr.d)) %>% slice_min(order_by = cdf.norm, n=n_top_clones)
            large.clones <- large.clones %>% filter(donor!=as.character(curr.d)) %>% bind_rows(new.large)
        }

    }
    return(large.clones)
}

get.top.clones <- function(d, clones, cdf_thresh, n_top_clones){
    curr.clones = clones %>% filter(donor==as.character(d))
    curr.clones = curr.clones %>% slice_min(order_by = cdf.norm, n=n_top_clones)
    cdf.clones = curr.clones %>% filter(cdf.norm<cdf_thresh)
    print('cdf thresh clones')
    print(dim(cdf.clones))

    if (dim(cdf.clones[1])>n_top_clones) {
        print("cdf too small. Taking top clones")
        print(n_top_clones)
        curr.clones = cdf.clones
    }
    return(curr.clones)
    }



plot.DE.RNA.pair <- function(integrated, de.results, a, b, outdir){
    try
    {
    plot1 <- VlnPlot(
      object = integrated,
      features = rownames(de.results)[1],
      pt.size = 0.1,
      idents = c(a,b)
    )
    plot2 <- FeaturePlot(
      object = integrated,
      features = rownames(de.results)[1],
      pt.size = 0.1
    )

    plot3 <- FeaturePlot(
      object = integrated,
      features = rownames(de.results)[2],
      pt.size = 0.1
    )

    plot1 | plot2 | plot3

    ggsave(file.path(outdir,paste0("clones_",a,"__",b, ".DE.GeneActivity.top2.png")))

    return(c(plot1, plot2, plot3))
    }
    return
}


## Run TF DE summary on filtered clone set
filter.clone.size <- function(clone.sizes, cdf.thresh){
    large.clones <- clone.sizes %>% filter(cdf.norm<cdf.thresh)
    return(large.clones)
}


setup.donor.clones <- function(se, d, large.clones, outdir){
        curr.outdir <- file.path(outdir, paste0("donor", d, "_TF"))
        dir.create(curr.outdir)
        donor.large.clones <- large.clones %>% filter(donor==d)
        clones.filt.ids <- donor.large.clones$lineage
        se.filt <- subset(se, subset = (donor==d) & (lineage %in% donor.large.clones$lineage))
        pairs = combn(clones.filt.ids,2)
        return(se.filt, pairs, curr.outdir)
}



#run.tf.large.clones <- function(se, outdir, n_donors, cdf.thresh)
summarize.tf.large.clones <- function(se, outdir, d, cdf.thresh, all.names.sig=c()){
    DefaultAssay(se) <- "chromvar"
    print('donor')
    print(d)
    list[se.filt, pairs, curr.outdir] =  setup.donor.clones(se, d, large.clones, outdir)

    names.sig <- c()
    for (i in 1:dim(pairs)[2]){
        print(pairs[,i])
            print(paste("clones", i))
            a = pairs[1,i]
            b = pairs[2,i]
            curr.tf.da <- read.csv(file.path(curr.outdir,
                                             paste0("clones_",a,"__", b,".DE.TF.csv"))) %>%
                          arrange(p_val)
            names.sig <- c(names.sig, head(curr.tf.da$X, n_top_genes))
            all.names.sig <- c(all.names.sig, head(curr.tf.da$X, n_top_genes))
            }
    names.sig <- unique(names.sig)
    names.sig


    vln <- VlnPlot(se.filt,  features=names.sig)
    dot <- DotPlot(se.filt, features = names.sig) + RotatedAxis()
    ggsave(plot=dot,
           file=file.path(curr.outdir, ("dot.top.png")))
    ggsave(plot=vln,
           file=file.path(curr.outdir, ("violin.top.png")))

    # split by a vector
    pdf((file.path(curr.outdir, "heatmap.top.pdf")), width=8,height=8)
    #pdf(qq("heatmap.pdf"), width = 8, height = 8)
    heat <- ComplexHeatmap::Heatmap(GetAssayData(se.filt)[names.sig,],
            name = paste0("donor",d),
            column_split = se.filt[[]]["lineage"],
            show_column_names = FALSE, use_raster=TRUE
           )
    ComplexHeatmap::draw(heat)
    dev.off()

    return(all.names.sig )
}





get.pwm <- function(se, genome, out_f=""){
  # Created by: ras
  # Created on: 10/26/21
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(future)

  # plan("multiprocess", workers = 16)
  # #options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
  # options(future.globals.maxSize = 8000 * 1024^2)
  # extract position frequency matrices for the motifs
  pwm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
  )

  fa.file <- Rsamtools::FaFile(genome, index=sprintf("%s.fai", genome)) #,

  DefaultAssay(se) <- "ATAC"
  chrom.assay <- se[["ATAC"]]
  # add motif information
  chrom.assay <- AddMotifs(chrom.assay, genome = fa.file, pfm = pwm)
  se = SetAssayData(se, slot="motifs", Motifs(chrom.assay))
  se = RunChromVAR(se, genome=fa.file)

  ## Add motif chromvarnames
  DefaultAssay(se) <- "chromvar"
  chrom.var.names <- GetAssayData(se)
  DefaultAssay(se) <- "ATAC"
  motifs <- Motifs(se)
  row.names(chrom.var.names) <- sapply(row.names(chrom.var.names), function(x) {Motifs(se)@motif.names[[x]]})

  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  se[['chromvarnames']] <- CreateAssayObject(data = chrom.var.names)
  DefaultAssay(se) <- "chromvarnames"


  if (out_f != ""){
    writeRDS(file=out_f, obj=se)
  }
  return(se)
}

