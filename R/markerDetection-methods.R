
# makeMarkerDetectionMatrix ----

#' Identify Markers and Signatures Present in Individual Samples
#'
#' Create matrices of detection events for individual gene features or combination thereof (i.e., signatures).
#'
#' @rdname makeDetectionMatrices
#'
#' @details
#' The \code{makeMarkerDetectionMatrix} function declares a feature as detected if it is detected above a given threshold in a specific assay (e.g., count or UMI matrix).
#'
#' The \code{makeSignatureDetectionMatrix} function declares a signature (composed of one or more gene features) as detected if all the associated features are detected.
#'
#' The \code{makeMarkerProportionMatrix} function computes the proportion of samples positive for each marker in predefined clusters.
#'
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param markers A character vector, subset of \code{rownames(se)}.
#' @param threshold Value \emph{above which} the marker is considered detected.
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#'
#' @return
#' \describe{
#' \item{\code{makeMarkerDetectionMatrix}}{A logical matrix indicating the presence of each marker in each sample.}
#' \item{\code{makeSignatureDetectionMatrix}}{A logical matrix indicating the presence of each signature in each sample.}
#' }
#'
#' @export
#' @importFrom SummarizedExperiment assay
#'
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' nsamples <- 100
#' u <- matrix(rpois(20000, 2), ncol=nsamples)
#' rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
#' se <- SummarizedExperiment(assays=list(counts=u))
#'
#' gsc <- GeneSetCollection(list(
#'   GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
#'   GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
#' ))
#'
#' # Example usage ----
#' markerMatrix <- makeMarkerDetectionMatrix(se, unique(unlist(geneIds(gsc))))
#' signatureMatrix <- makeSignatureDetectionMatrix(markerMatrix, gsc)
makeMarkerDetectionMatrix <- function(
    se, markers, threshold=0, assay.type="counts"
) {
    if (any(duplicated(markers))) {
        warning("Dropping duplicated markers values")
        markers <- unique(markers)
    }
    # Subset the requested assay to the markers of interest
    x <- assay(se, assay.type)[markers, , drop=TRUE]
    markerDetectionMatrix <- (x > threshold)
    markerDetectionMatrix
}

# makeSignatureDetectionMatrix ----

#' @rdname makeDetectionMatrices
#'
#' @param matrix A logical matrix indicating the presence of each marker in each sample.
#' See \code{\link{makeMarkerDetectionMatrix}}
#' @param object A collection of signatures inheriting from "\code{\link{GeneSetCollection}}" or "\code{\link{tbl_geneset}}".
#'
#' @export
#'
#' @importFrom S4Vectors FilterRules evalSeparately
makeSignatureDetectionMatrix <- function(
    matrix, object
) {
    filterExpressions <- makeFilterExpression(object)
    fr <- FilterRules(filterExpressions)
    es <- evalSeparately(fr, as.data.frame(t(matrix)))
    es
}

# makeMarkerProportionMatrix ----

#' @rdname makeDetectionMatrices
#'
#' @param cluster.col Name of a column in \code{colData(se)} that contains
#' a factor indicating cluster membership for each column (i.e. sample) in \code{se}.
#'
#' @export
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom Matrix rowSums
makeMarkerProportionMatrix <- function(
    se, cluster.col, assay.type="counts", threshold=0
) {
    stopifnot(!missing(cluster.col))
    stopifnot(is.factor(colData(se)[, cluster.col, drop=TRUE])) # keep as-is to raise an error referring to variables known to the user
    clusterData <- colData(se)[, cluster.col, drop=TRUE]

    # Compute the proportion of each cluster positive for each marker
    markerDetectionMatrix <- makeMarkerDetectionMatrix(se, rownames(se), threshold, assay.type)

    clusterNames <- levels(clusterData)
    numberCellsInCluster <- table(clusterData)

    proportionPositiveByCluster <- matrix(
        data=NA_real_,
        nrow=nrow(markerDetectionMatrix),
        ncol=length(clusterNames),
        dimnames=list(feature=rownames(markerDetectionMatrix), cluster=clusterNames))
    x <- assay(se, assay.type)
    for (clusterName in clusterNames) {
        clusterSamples <- which(colData(se)[, cluster.col] == clusterName)
        nDetected <- Matrix::rowSums(x[, clusterSamples] > threshold)
        proportionPositiveByCluster[, clusterName] <- nDetected / length(clusterSamples)
    }

    proportionPositiveByCluster
}
