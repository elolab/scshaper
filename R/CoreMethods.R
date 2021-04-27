#' @title Run \code{scShaper}.
#'
#' @description
#' This function wraps all steps in scShaper from dimensionality reduction to LOESS smoothing.
#'
#' @param object an object of \code{SingleCellExperiment} class
#' @param dim.red.method a character denoting which dimensionality reduction method to use.
#' Default is \code{tsne}. Other options are PCA (\code{pca}) 
#' and no dimensionality reduction (\code{none}).
#' @param num.pcs a positive integer denoting how many principal components (PCs) are extracted. 
#' Default is \code{50}.
#' @param num.tsne a positive integer denoting how many t-SNE dimensions to extract. 
#' Default is \code{3}.
#' @param tsne.perplexity a numeric denoting the perplexity of t-SNE. Default is \code{30}.
#' @param k.range a numeric vector specifying the different cluster numbers, 
#' for which the discrete pseudotime is estimated using k-means and 
#' the special case of Kruskal's algorithm.
#' Default is \code{2:100}. 
#' The first parameter to consider increasing if unsatisfactory results are obtained.
#' @param iter.max a positive integer denoting how iterations to perform 
#' at maximum in k-means. Default is \code{1000}.
#' @param nstart a positive integer denoting how many initializations to perform 
#' in k-means. Default is \code{1}. 
#' The third parameter to consider increasing if unsatisfactory results are obtained.
#' @param pca.rank a positive integer denoting how many principal components to consider 
#' when generating the ensemble pseudotime. Default is \code{2}. 
#' @param span a numeric denoting the degree of smoothing in LOESS. 
#' Default is \code{0.1}. The second parameter to consider increasing if unsatisfactory results are obtained.
#' @name RunscShaper
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords Seurat graph clustering PCA feature selection
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom methods is
#'
#' @examples
#' library(SingleCellExperiment)
#' library(scShaper)
#' sce <- SingleCellExperiment(assays = list(counts=prosstt1$rna.raw,logcounts = prosstt1$rna.normalized))
#' sce <- RunscShaper(sce)

RunscShaper.SingleCellExperiment <- function(object,dim.red.method,num.pcs,num.tsne,tsne.perplexity,k.range,iter.max,nstart,pca.rank,span) {
  
  # Prepare scShaper
  object <- PreparescShaper(object)
  
  # Run dimensionality reduction
  object <- RunDimRed(object,dim.red.method,num.pcs,num.tsne,tsne.perplexity)
  
  message("1. Run dimensionality reduction. READY")
  
  # Run k-means clustering
  object <- RunKMeans(object,k.range,iter.max,nstart)

  message("2. Run clustering. READY")

  # Run special case of Kruskal's algorithm
  object <- RunSHP(object)
  
  message("3. Run the special case of Kruskal's algorithm. READY")
  
  # 4. Run a wrapper function for part of the analysis that finds the largest subset
  # of the discrete pseudotimes that are mutually correlated.
  object <- RunPCAAnalysis(object,pca.rank=pca.rank)
  
  # To generate the final continuous pseudotime, run epsilon SVR to smooth the discrete pseudotime
  message("4. Generate average pseudotime based on the largest mutually correlating subset of the different discrete pseudotimes. READY")

  object <- RunLoessSmoothing(object,span)
  
  message("5. Run LOESS based smoothing. READY")
  return(object)
}

#' @rdname RunscShaper
#' @aliases RunscShaper
setMethod("RunscShaper", signature(object = "SingleCellExperiment"),
          RunscShaper.SingleCellExperiment)



#' @title Prepare \code{SingleCellExperiment} object for \code{scShaper} analysis
#'
#' @description
#' This function prepares the \code{SingleCellExperiment} object for
#' \code{scShaper} analysis. The only required input is an object of class
#' \code{SingleCellExperiment} with data stored at least in the \code{logcounts} and \code{counts} slots.
#'
#' @param object an object of \code{SingleCellExperiment} class
#' @param sparse.data a logical denoting if the data should
#' be converted into sparse format to save memory. Default is \code{TRUE}.
#'
#' @name PreparescShaper
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords prepare scshaper raw normalized data
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom methods is
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts=prosstt1$rna.raw,logcounts = prosstt1$rna.normalized))
#' sce <- PreparescShaper(sce)
#'
PreparescShaper.SingleCellExperiment <- function(object,sparse.data) {
  
  # Check that there are data in `logcounts` slot
  if (!("logcounts" %in% assayNames(object))) {
    stop(paste("`Error: `logcounts` slot is missing from your ",
               "SingleCellExperiment object. This can be any kind of ",
               "normalized data matrix. Set it by executing ",
               "logcounts(object) <- norm_data",sep = ""))
    return(object)
  }
  
  # Check that there are data in `counts` slot
  if (!("counts" %in% assayNames(object))) {
    stop(paste("`Error: `counts` slot is missing from your ",
               "SingleCellExperiment object. This can be any kind of ",
               "normalized data matrix. Set it by executing ",
               "counts(object) <- raw_data",sep = ""))
    return(object)
  }
  
  
  # Remove duplicate features
  if (sum(duplicated(rownames(object))) != 0) {
    features_before <- length(rownames(object))
    object <- object[!duplicated(rownames(object)), ]
    features_after <- length(rownames(object))
    message(paste("data in SingleCellExperiment object contained duplicate ",
                  " features. ", features_before - features_after,
                  "/", features_before, " were filtered out."))
  }
  
  # Convert the data in `logcounts` slot into object of `dgCMatrix` class.
  if (is(logcounts(object), "matrix") & sparse.data) {
    logcounts(object) <- Matrix(logcounts(object),sparse = TRUE)
    message(paste("Converting object of `matrix` class into `dgCMatrix`.",sep=""))
  }
  else if (is(logcounts(object), "data.frame") & sparse.data) {
    logcounts(object) <- Matrix(as.matrix(logcounts(object)),sparse = TRUE)
    message(paste("Converting object of `data.frame` class into `dgCMatrix`.",
                  " Please note that scShaper has been designed to work with ",
                  "sparse data, i.e. data with ",
                  "a high proportion of zero values!",sep = ""))
  }
  else if (is(logcounts(object), "dgCMatrix")) {
    message("Data in `logcounts` slot already of `dgCMatrix` class...")
  }
  else {
    stop("Error: Data in `logcounts` slot is not of `matrix`, `data.frame` ",
         "or `dgCMatrix` class.")
    return(object)
  }
  
  
  # Convert the data in `counts` slot into object of `dgCMatrix` class.
  if (is(counts(object), "matrix") & sparse.data) {
    counts(object) <- Matrix(counts(object),sparse = TRUE)
    message(paste("Converting object of `matrix` class into `dgCMatrix`.",sep=""))
  }
  else if (is(counts(object), "data.frame") & sparse.data) {
    counts(object) <- Matrix(as.matrix(counts(object)),sparse = TRUE)
    message(paste("Converting object of `data.frame` class into `dgCMatrix`.",sep = ""))
  }
  else if (is(counts(object), "dgCMatrix")) {
    message("Data in `counts` slot already of `dgCMatrix` class...")
  }
  else {
    stop("Error: Data in `counts` slot is not of `matrix`, `data.frame` ",
         "or `dgCMatrix` class.")
    return(object)
  }
  
  
  # Filter genes that are not expressed in any of the cells
  genes_before_filtering <- nrow(object)
  non_expressing_genes <- rownames(object)[rowSums(logcounts(object)) != 0]
  object <- object[non_expressing_genes,]
  genes_after_filtering <- nrow(object)
  message(paste(genes_after_filtering,"/",genes_before_filtering,
                " genes remain after filtering genes with only zero values.",
                sep = ""))
  
  # Create a place into `metadata`` slot for the data from scShaper
  metadata(object)$scshaper <- list()
  
  return(object)
}

#' @rdname PreparescShaper
#' @aliases PreparescShaper
setMethod("PreparescShaper", signature(object = "SingleCellExperiment"),
          PreparescShaper.SingleCellExperiment)



#' @title Run dimensionality reduction.
#'
#' @description
#' This function performs dimensionality reduction using just PCA or PCA and t-SNE (default).
#' @param object an object of \code{SingleCellExperiment} class
#' @param dim.red.method a character denoting which dimensionality reduction method to use.
#' Default is \code{tsne}. Other options are PCA (\code{pca}) 
#' and no dimensionality reduction (\code{none}).
#' @param num.pcs a positive integer denoting how many principal components (PCs) are extracted. 
#' Default is \code{50}.
#' @param num.tsne a positive integer denoting how many t-SNE dimensions to extract. 
#' Default is \code{3}.
#' @param tsne.perplexity a numeric denoting the perplexity of t-SNE. Default is \code{30}.
#' @name RunDimRed
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords Seurat graph clustering PCA feature selection
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<- reducedDim
#' @importFrom methods is
#' @importFrom Rtsne Rtsne
#'
#' @examples
#' library(SingleCellExperiment)
#' library(scShaper)
#' sce <- SingleCellExperiment(assays = list(counts=prosstt1$rna.raw,logcounts = prosstt1$rna.normalized))
#' sce <- PreparescShaper(sce)
#' sce <- RunDimRed(sce)

RunDimRed.SingleCellExperiment <- function(object,dim.red.method,num.pcs,num.tsne,tsne.perplexity) {
  
  
  X <- Matrix::t(logcounts(object))
  
  if (dim.red.method=="tsne")
  {
    if (nrow(X) < 90)
    {
      tsne.perplexity <- 20
    }
    else if (nrow(X) < 60)
    {
      tsne.perplexity <- 15
    }
    rtsne_out <- Rtsne(as.matrix(X),dims=num.tsne,
                       initial_dims=num.pcs,
                       perplexity=tsne.perplexity)
    X <- rtsne_out$Y
    reducedDim(object,type = "tsne") <- X
  }
  
  else if (dim.red.method=="pca")
  {
    pca_out <- prcomp(X,scale. = TRUE,rank. = num.pcs)
    X <- pca_out$x[,1:num.pcs]
    reducedDim(object,type = "pca") <- X
  }
  else {
    message("No dimensionality reduction performed!")
  }
  
  return(object)
}

#' @rdname RunDimRed
#' @aliases RunDimRed
setMethod("RunDimRed", signature(object = "SingleCellExperiment"),
          RunDimRed.SingleCellExperiment)





#' @title Run k-means clustering
#'
#' @description
#' This function runs k-means clustering for a range of k values.
#' @param object an object of \code{SingleCellExperiment} class
#' @param k.range a numeric vector specifying the different cluster numbers, 
#' for which the discrete pseudotime is estimated using k-means and 
#' the special case of Kruskal's algorithm.
#' Default is \code{2:100}. 
#' @param iter.max a positive integer denoting how iterations to perform 
#' at maximum in k-means. Default is \code{1000}.
#' @param nstart a positive integer denoting how many initializations to perform 
#' in k-means. Default is \code{1}. 
#' The third parameter to consider increasing if unsatisfactory results are obtained.
#' @name RunKMeans
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords Seurat graph clustering PCA feature selection
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<- reducedDim reducedDimNames
#' @importFrom methods is
#' @importFrom Rtsne Rtsne
#' @importFrom stats kmeans
#'
#' @examples
#' library(SingleCellExperiment)
#' library(scShaper)
#' sce <- SingleCellExperiment(assays = list(counts=prosstt1$rna.raw,logcounts = prosstt1$rna.normalized))
#' sce <- PreparescShaper(sce)
#' sce <- RunDimRed(sce)
#' sce <- RunKMeans(sce)

RunKMeans.SingleCellExperiment <- function(object,k.range,iter.max,nstart) {
  
  if (identical(reducedDimNames(object), character(0)))
  {
    message("No dimensionality reduction available. Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- Matrix::t(X)
  } else
  {
    X <- reducedDim(object)
  }
  
  if (max(k.range) >= nrow(X))
  {
    k.range <- seq(min(k.range),nrow(X)-1)
  }
  
  kmeans_out_list <- list()
  for (k in k.range)
  {
    kmeans_out <- kmeans(X,centers = k,iter.max = iter.max,nstart = nstart)
    kmeans_out_list[[as.character(k)]] <- kmeans_out$cluster
  }
  metadata(object)$scshaper$clustering <- kmeans_out_list
  
  return(object)
}

#' @rdname RunKMeans
#' @aliases RunKMeans
setMethod("RunKMeans", signature(object = "SingleCellExperiment"),
          RunKMeans.SingleCellExperiment)



#' @title Run a special case of kruskal's algorithm to find an optimal path through each clustering
#'
#' @description
#' Run a special case of kruskal's algorithm to find an optimal path through each clustering:
#' shortest Hamiltonian path (SHP)
#' @param object an object of \code{SingleCellExperiment} class
#' @name RunSHP
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords Seurat graph clustering PCA feature selection
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<- reducedDim
#' @importFrom methods is
#' @importFrom plyr mapvalues
#'
#' @useDynLib scShaper
#'
#' @examples
#' library(SingleCellExperiment)
#' library(scShaper)
#' sce <- SingleCellExperiment(assays = list(counts=prosstt1$rna.raw,logcounts = prosstt1$rna.normalized))
#' sce <- PreparescShaper(sce)
#' sce <- RunDimRed(sce)
#' sce <- RunKMeans(sce)
#' sce <- RunSHP(sce)

RunSHP.SingleCellExperiment <- function(object) {
  
  if (identical(reducedDimNames(object), character(0)))
  {
    message("No dimensionality reduction available. Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- Matrix::t(X)
  } else
  {
    X <- reducedDim(object)
  }
  
  # Extract clustering
  clustering <- metadata(object)$scshaper$clustering
  
  msp_out_list <- list()
  for (k in names(clustering))
  {
    clustering_k <- clustering[[k]]
    # Sort clustering using MSP algorithm
    cluster_centers <- do.call(rbind,lapply(split(1:nrow(X),clustering_k),function(x) apply(X[x,,drop=FALSE],2,mean)))
    dm <- as.matrix(dist(cluster_centers))
    clusters_sorted <- specialKruskal(dm)
    clustering_sorted <- mapvalues(clustering_k,from = as.numeric(names(table(clustering_k))),order(as.numeric(unlist(strsplit(clusters_sorted,"_")))))
    msp_out_list[[k]] <- clustering_sorted
  }
  metadata(object)$scshaper$discrete.pseudotime <- msp_out_list
  
  return(object)
}

#' @rdname RunSHP
#' @aliases RunSHP
setMethod("RunSHP", signature(object = "SingleCellExperiment"),
          RunSHP.SingleCellExperiment)







#' @title Run PCA
#'
#' @description
#' This function performs PCA to analyze linear dependencies between the discrete pseudotimes.
#' @param object an object of \code{SingleCellExperiment} class
#' @param pca.rank a positive integer denoting how many principal components to consider 
#' when generating the ensemble pseudotime. Default is \code{2}. 
#' @name RunPCAAnalysis
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords Seurat graph clustering PCA feature selection
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<- reducedDim
#' @importFrom methods is
#' @importFrom plyr mapvalues
#' @importFrom stats prcomp
#'
#' @examples
#' library(SingleCellExperiment)
#' library(scShaper)
#' sce <- SingleCellExperiment(assays = list(counts=prosstt1$rna.raw,logcounts = prosstt1$rna.normalized))
#' sce <- PreparescShaper(sce)
#' sce <- RunDimRed(sce)
#' sce <- RunKMeans(sce)
#' sce <- RunSHP(sce)
#' sce <- RunPCAAnalysis(sce)

RunPCAAnalysis.SingleCellExperiment <- function(object,pca.rank) {
  
  
  dp <- do.call(cbind,object@metadata$scshaper$discrete.pseudotime)

  
  discrete.pseudotime <- metadata(object)$scshaper$discrete.pseudotime


  # Combine all discrete pseudotimes to column matrix
  dptm <- do.call(cbind,discrete.pseudotime)
  # Min-max scale to range [0,1]
  dptm <- apply(dptm,2,function(x) (x-min(x))/(max(x)-min(x)))

  # Perform PCA
  pca_out <- prcomp(dp,scale. = TRUE,rank. = pca.rank)
  
  # Select the discrete pseudotimes that contribute most to the 1st PC
  clusterings_selected <- apply(abs(pca_out$rotation[,1:pca.rank]),1,which.max) == 1
  dpt <- dptm[,clusterings_selected]
  
  # Flip half of the pseudotimes
  sign_selected <- sign(pca_out$rotation[,1])[clusterings_selected]
  
  if (sum(sign_selected==-1)!=0)
  {
    dpt_to_flip <- dpt[,sign_selected==-1,drop=FALSE]
    dp_flipped <- apply(dpt_to_flip,2,function(x) 1-x)
    dpt[,sign_selected==-1] <- dp_flipped
  }
  
  
  # Calculate average of the discrete pseudotimes
  metadata(object)$scshaper$average.discrete.pseudotime <- apply(dpt,1,mean)
  
  return(object)
}

#' @rdname RunPCAAnalysis
#' @aliases RunPCAAnalysis
setMethod("RunPCAAnalysis", signature(object = "SingleCellExperiment"),
          RunPCAAnalysis.SingleCellExperiment)







#' @title Run LOESS smoothing
#'
#' @description
#' This function performs LOESS to smoothen the crude initial average discrete pseudotime.
#' @param object an object of \code{SingleCellExperiment} class
#' @param span a numeric denoting the degree of smoothing in LOESS. 
#' Default is \code{0.1}. The second parameter to consider increasing 
#' if unsatisfactory results are obtained.
#' @name RunLoessSmoothing
#'
#' @return an object of \code{SingleCellExperiment} class
#'
#' @keywords Seurat graph clustering PCA feature selection
#'
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Matrix rowSums Matrix
#' @importFrom SingleCellExperiment logcounts logcounts<- reducedDim
#' @importFrom methods is
#' @importFrom stats loess
#' @examples
#' library(SingleCellExperiment)
#' library(scShaper)
#' sce <- SingleCellExperiment(assays = list(counts=prosstt1$rna.raw,logcounts = prosstt1$rna.normalized))
#' sce <- PreparescShaper(sce)
#' sce <- RunDimRed(sce)
#' sce <- RunKMeans(sce)
#' sce <- RunSHP(sce)
#' sce <- RunPCAAnalysis(sce)
#' sce <- RunLoessSmoothing(sce)

RunLoessSmoothing.SingleCellExperiment <- function(object,span) {
  
  if (identical(reducedDimNames(object), character(0)))
  {
    message("No dimensionality reduction available. Using data in the 'logcounts' slot.")
    X <- logcounts(object)
    X <- Matrix::t(X)
  } else
  {
    X <- reducedDim(object)
  }
  
  # Extract clustering
  discrete_pseudotime <- metadata(object)$scshaper$average.discrete.pseudotime
  
  df <- as.data.frame(cbind(rank(discrete_pseudotime),discrete_pseudotime))
  colnames(df) <- paste0("V",1:2)
  lo <- loess(V2 ~ V1,df,span = span,degree = 2)
  metadata(object)$scshaper$continuous.pseudotime <- predict(lo) 
  
  
  
  return(object)
}

#' @rdname RunLoessSmoothing
#' @aliases RunLoessSmoothing
setMethod("RunLoessSmoothing", signature(object = "SingleCellExperiment"),
          RunLoessSmoothing.SingleCellExperiment)



