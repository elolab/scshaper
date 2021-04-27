#' @export
setGeneric("PreparescShaper",signature = "object",
           function(object,sparse.data=TRUE) {
             standardGeneric("PreparescShaper")
           })

#' @export
setGeneric("RunscShaper",signature = "object",
           function(object,dim.red.method="tsne",num.pcs=50,num.tsne=3,tsne.perplexity=30,k.range=seq(2,100),iter.max = 1000,nstart = 1,pca.rank=2,span=0.1) {
             standardGeneric("RunscShaper")
           })

#' @export
setGeneric("RunDimRed",signature = "object",
           function(object,dim.red.method="tsne",num.pcs=50,num.tsne=3,tsne.perplexity=30) {
             standardGeneric("RunDimRed")
           })

#' @export
setGeneric("RunKMeans",signature = "object",
           function(object,k.range=seq(10,100),iter.max=1000,nstart=1) {
             standardGeneric("RunKMeans")
           })

#' @export
setGeneric("RunSHP",signature = "object",
           function(object) {
             standardGeneric("RunSHP")
           })

#' @export
setGeneric("RunPCAAnalysis",signature = "object",
           function(object,pca.rank=5) {
             standardGeneric("RunPCAAnalysis")
           })


#' @export
setGeneric("RunLoessSmoothing",signature = "object",
           function(object,span=0.1) {
             standardGeneric("RunLoessSmoothing")
           })

