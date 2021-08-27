library(Matrix)
library(SummarizedExperiment)
#install.packages("tidyverse")
#install.packages("umap")
library(tidyverse)
library(umap)
library(edgeR)
library(FNN)
library(matrixStats)
library(igraph)
library(BuenColors)
library(GenomicRanges)
set.seed(1)
"%ni%" <- Negate("%in%")


#Binarize Sparse Matrix
binarizeMat <- function(mat){
  mat@x[mat@x > 0] <- 1
  mat
}

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL, log = TRUE){
  
  set.seed(1)
  
  #TF IDF LSI adapted from flyATAC
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
    mat <- mat[idx,] 
  }else{
    idx <- which(Matrix::rowSums(mat) > 0)
    mat <- mat[idx,]
  }
  
  #Calc RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm)
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
  if(log){
    freqs@x = log1p(freqs@x * 100000) # equivalent to adding a small pseudocount, but without making matrix dense
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
    
  } else{
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  }
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  #Return Object
  out <- list(
    matSVD = matSVD, 
    rowSm = rowSm, 
    colSm = colSm, 
    idx = idx, 
    svd = svd, 
    binarize = binarize, 
    nComponents = nComponents,
    date = Sys.Date(),
    seed = 1)
  
  out
  
}

projectLSI <- function(mat, lsi, log = TRUE){   
  
  #Get Same Features
  mat <- mat[lsi$idx,]
  if(lsi$binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1       
  }
  
  #Calc TF IDF
  rowsToZero <- which(lsi$rowSm == 0)
  setToZero <- which((mat@i + 1) %in% rowsToZero)
  if(length(setToZero) > 0){
    mat@x[setToZero] <- 0
  }
  
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
  if(log){
    freqs@x = log1p(freqs@x * 100000) # equivalent to adding a small pseudocount, but without making matrix dense
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
    
  } else{
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  }
  if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
    tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
  }
  
  #Calc V
  V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
  diag(svdDiag) <- lsi$svd$d
  matSVD <- t(svdDiag %*% t(V))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  return(matSVD)
  
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

#Seurat SNN
seuratSNN <- function(matSVD, dims.use = 1:50, print.output = TRUE, ...){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  obj <- Seurat::CreateSeuratObject(mat, project='scATAC', min.cells=0, min.genes=0)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
  obj <- Seurat::FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, print.output = print.output, ...)
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  paste0("mc",match(clust, unique(clust)))
}

louvainIgraphClusters <- function(matSVD, k){
  set.seed(1)
  knn <- FNN::get.knn(matSVD, algo="kd_tree", k = k)[["nn.index"]]
  igraphObj <- igraph::graph_from_adjacency_matrix(igraph::get.adjacency(igraph::graph.edgelist(data.matrix(reshape2::melt(knn)[,c("Var1", "value")]), directed=FALSE)), mode = "undirected")
  clusters <- membership(cluster_louvain(igraphObj))
  return(paste0("mc", as.character(clusters)))
}

build_means_for_clusters <- function(LSI_mat, clusters){
  sapply(sort(unique(clusters)), function(cluster_i){
    colMeans(LSI_mat[clusters == cluster_i,])
  })
}

# A: A matrix of clusters x PCs
# B: A matrix of cells x PCs
classify_from_reference <- function(A, B){
  
  # Get the pair-wise euclidean distance most-alike
  euklDist <- t(sqrt(apply(array(apply(B,1,function(x){(x-t(A))^2}),c(ncol(A),nrow(A),nrow(B))),2:3,sum)))
  colnames(euklDist) <- rownames(A)
  colnames(euklDist)[max.col(-1*euklDist, 'first')]  -> vec
  return(vec)
}
