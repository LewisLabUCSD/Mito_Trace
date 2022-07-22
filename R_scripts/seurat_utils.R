library(Signac)
library(Seurat)
library(dplyr)
library("Nebulosa")


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


embed.atac <- function(se, outdir, lsi_start_comp=2, lsi_end_comp=50,reduction='lsi',
                       neighbor_dim=30,return.depth=T){
# Binarize and run LSI
    se <- FindTopFeatures(se, min.cutoff = 20)
    se <- BinarizeCounts(se)
    se <- RunTFIDF(se)
    se <- RunSVD(se)
    se <- RunUMAP(se, dims = lsi_start_comp:lsi_end_comp, reduction = reduction)
    #dimP <- DimPlot(se, group.by = "proj", pt.size = 0.1)
    pDepthCorr <- DepthCor(se, reduction=reduction)
    #ggsave(file.path(outdir,"integrated.depthCor.png"), plot=pDepthCorr, dpi=300)
    #pDepthCorr
        #integrated <- RunUMAP(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
    se <- FindNeighbors(object = se, reduction = 'lsi', dims = lsi_start_comp:lsi_end_comp)
    se <- FindClusters(object = se, verbose = FALSE, algorithm = 3)
    if (return.depth){
      return(c(se, pDepthCorr))
    }else{return(se)}
}


featplot <- function(name.sig, se, curr.outdir, feat.names=NULL, use.density=FALSE){
  if (!is.null(feat.names)){
    name <- feat.names[name.sig,]
  }else{
    name <- name.sig
  }
  if (name.sig %in% rownames(se)){
    feat <- FeaturePlot(se,  features=name.sig) + ggtitle(name)
    ggsave(plot=feat,
           file=file.path(curr.outdir, paste0(name.sig,".embedFeat.top.png")))
  }else{
    print(paste0("Feature no in object: ", name.sig))
  }
  if(use.density){
    plot_density(se, name.sig) + ggtitle(name)
    ggsave(plot=feat,
           file=file.path(curr.outdir, paste0("density.", name.sig,".embedFeat.top.png")))
  }
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
            new.large <- init.large.clones %>% dplyr::filter(donor==as.character(curr.d)) %>% slice_min(order_by = cdf.norm, n=n_top_clones)
            large.clones <- large.clones %>% dplyr::filter(donor!=as.character(curr.d)) %>% bind_rows(new.large)
        }

    }
    return(large.clones)
}

get.top.clones <- function(d, clones, cdf_thresh, n_top_clones){
    curr.clones = clones %>% dplyr::filter(donor==as.character(d))
    curr.clones = curr.clones %>% dplyr::slice_min(order_by = cdf.norm, n=n_top_clones)
    cdf.clones = curr.clones %>% dplyr::filter(cdf.norm<cdf_thresh)
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


de.plots <- function(se.filt, names.sig, curr.outdir, curr.name="", max.size=10, to.heat=T,to.vln=T,
                     feature.names=NULL){
    if (length(names.sig) > max.size){
        names.sig <- names.sig[1:max.size]
    }
    if (is.null(feature.names)){
      names <- names.sig
    }else{
      names <- feature.names[1:length(names.sig)]
    }
    
    dot <- DotPlot(se.filt, features = names.sig) + RotatedAxis() + ggtitle(curr.name) + scale_x_discrete(labels=names) 
    
    feat <- FeaturePlot(se.filt,  features=names.sig)
    if (to.vln){
    vln <- VlnPlot(se.filt,  features=names.sig, pt.size = 0)
    }
    # split by a vector
    if (to.heat){
    pdf(file.path(curr.outdir, paste0(curr.name, ".heatmap.top.pdf")), width=8,height=8)

    heat <- ComplexHeatmap::Heatmap(as.matrix(GetAssayData(se.filt)[names.sig,]),
            name = curr.name,
            show_column_names = FALSE, use_raster=TRUE
           )
    ComplexHeatmap::draw(heat)
    dev.off()
    }
    ggsave(plot=feat,
           file=file.path(curr.outdir, paste0(curr.name,".embedFeat.top.png")))
    ggsave(plot=dot,
           file=file.path(curr.outdir, paste0(curr.name, ".dot.top.png")))
    if (to.vln){
    ggsave(plot=vln,
           file=file.path(curr.outdir, paste0(curr.name, ".violin.top.png")))
           ggsave(plot=vln,
                  file=file.path(curr.outdir, paste0(curr.name, ".violin.top.pdf")))
           }
    ## pdfs
    ggsave(plot=dot,
           file=file.path(curr.outdir, paste0(curr.name, ".dot.top.pdf")))
    return(dot)

}


find.markers.and.plot <- function(se, id1, id2, curr.outdir, curr.name, min.pct, logfcthresh,
                                  p.thresh=0.1,
                                  latent.vars=NULL, test.use="wilcox", assay="RNA" ,
                                  to.plot=T, fc.name=NULL, mean.fxn=NULL){
    print('se filt before')
    se.filt <- subset(se, idents = c(id1,id2))
    print('se filt')
    print(se.filt)
    response <- FindMarkers(se,
                            ident.1 = id1,
                            ident.2 = id2,
                            verbose = T,
                            only.pos = FALSE,
                            mean.fxn = mean.fxn,
                            logfc.threshold = logfcthresh,
                            min.pct = min.pct, test.use=test.use,
                            latent.vars=latent.vars,
                            fc.name = fc.name)
              
    ncells <- data.frame(table(Idents(se.filt)))
    write.csv(ncells, file=file.path(curr.outdir, paste0(curr.name,".counts.csv")), quote=F)
    response$p_val_adj_BH <- stats::p.adjust(response$p_val, method = "BH", n = length(response$p_val))
    response <- response %>% dplyr::arrange(p_val_adj_BH)
    
    curr.sig <- response %>% dplyr::filter(p_val_adj_BH<p.thresh)
    print("curr sig")
    print(head(curr.sig))
    #names.sig <- names(curr.sig)
    names.sig <- rownames(curr.sig)
    
    if (assay == "ATAC"){
      response.features <- ClosestFeature(se.filt, regions = rownames(response))
      response.features = response.features %>% dplyr::mutate(gene.id = paste(gene_name,type, sep="_"))
      response$gene.id = response.features$gene.id
    }
    
    #names.sig <- clean.de(response, se, n_top_genes, a=id1, b=id2, names.sig = c())
    #response <- curr.de[[1]]
    #names.sig <- curr.de[[2]]
    print('dim response')
    print(dim(response))
    if (!(dim(response)[1]==0)){
        write.csv(response,
                  file=file.path(curr.outdir, paste(curr.name,"DE.csv",sep="_")), quote=F)
      if (!(length(names.sig)==0)){
        # print('response plots')
        # print(head(response))
        # print(head(names.sig))
        if (to.plot){
          if (assay == "ATAC"){
            de.plots(se.filt, names.sig, curr.outdir, curr.name=curr.name, feature.names=response[names.sig,]$gene.id)
          }else{
            de.plots(se.filt, names.sig, curr.outdir, curr.name=curr.name)
          }
          try({
              print("gally")
              print(colnames(response))
              if (is.null(mean.fxn)){
                mean.fxn = "avg_log2FC"
              }
              gally <- GGally::ggpairs(response[,c("p_val", "p_val_adj", mean.fxn, "p_val_adj_BH")], aes(alpha = 0.4))
              ggsave(plot=gally, file=file.path(curr.outdir,
                                                paste(curr.name,".DE.pvalHist.png", sep="_")))
          })
        }
      }
    }
    return(curr.sig)
}

## Run TF DE summary on filtered clone set
filter.clone.size <- function(clone.sizes, cdf.thresh){
    large.clones <- clone.sizes %>% dplyr::filter(cdf.norm<cdf.thresh)
    return(large.clones)
}


setup.donor.clones <- function(se, d, large.clones, outdir){
        curr.outdir <- file.path(outdir, paste0("donor", d, "_TF"))
        dir.create(curr.outdir)
        donor.large.clones <- large.clones %>% dplyr::filter(donor==d)
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



# process.cells.meta <- function(se, cells_meta){
#   new_cells_meta = merge(cells_meta, se[[]], by=0, all=TRUE)
#   rownames(new_cells_meta) <- new_cells_meta$Row.names
#   head(new_cells_meta)
#   new_cells_meta$name <- paste0(new_cells_meta$donor, "_", new_cells_meta$lineage)
#   
#   donor_levels <- levels(factor(new_cells_meta$donor))
#   donor_levels[length(donor_levels) + 1] <- "None"
#   
#   condition_levels <- levels(factor(new_cells_meta$condition))
#   condition_levels[length(condition_levels) + 1] <- "None"
#   
#   lineage_levels <- levels(factor(new_cells_meta$lineage))
#   lineage_levels[length(lineage_levels) + 1] <- "None"
#   
#   new_cells_meta$name <- paste0(new_cells_meta$donor, "_", new_cells_meta$lineage)
#   name_levels <- levels(factor(new_cells_meta$name))
#   name_levels[length(name_levels) + 1] <- "None"
#   
#   new_cells_meta$donor <- factor(new_cells_meta$donor, levels=donor_levels)
#   new_cells_meta$condition <- factor(new_cells_meta$condition, levels=condition_levels)
#   new_cells_meta$lineage <- factor(new_cells_meta$lineage, levels=lineage_levels)
#   new_cells_meta$name <- factor(new_cells_meta$name, levels=name_levels)
#   
#   new_cells_meta$donor[is.na(new_cells_meta$donor)] <- "None"
#   new_cells_meta$condition[is.na(new_cells_meta$condition)] <- "None"
#   new_cells_meta$lineage[is.na(new_cells_meta$lineage)] <- "None"
#   new_cells_meta$name[(is.na(new_cells_meta$name) | grepl("NA", new_cells_meta$name))] <- "None"
#   return(new_cells_meta)
# }

add.cell.meta <- function(se, new_cells_meta){
  #se$ID <- apply(stringr::str_split((rownames(se[[]])), "_", simplify=T)[,1:2], 1, function(x) {paste0(x[[1]],  "_", x[[2]])})
  #se <- RenameCells(se, new.names=se$ID)
  se <- AddMetaData(se, new_cells_meta)
  se
  return(se)
}


get.clone.sizes <- function(se){
  counts <- se[[]] %>%
    group_by(condition, donor, lineage, name) %>%
    summarize(size=n()) %>% dplyr::filter(!name=="None")
  
  counts.norm <- counts %>% group_by(condition, donor) %>% mutate(total=sum(size)) %>% ungroup() %>% mutate(norm=size/total)
  
  counts.conds.norm = counts.norm[order(counts.norm[,"norm"], decreasing=T),]
  counts.conds.norm
  clone.sizes <- counts.conds.norm %>% 
    group_by(donor,lineage) %>% 
    summarize(norm.total=sum(norm), total=sum(size)) %>% 
    arrange(desc(norm.total)) %>% 
    mutate(name=factor(paste0(donor,"_", lineage)))
  
  clone.sizes <- clone.sizes %>% group_by(donor) %>% arrange(donor, desc(norm.total)) %>% 
    mutate(cdf.norm=cumsum(norm.total)/2, cdf =cumsum(total),
           index=1:n()) %>% 
    ungroup %>% arrange(donor, desc(norm.total)) 
  
  return(list("counts.norm"=counts.conds.norm, "clone.size"= clone.sizes))
}
