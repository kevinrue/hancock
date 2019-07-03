
# makeMarkerDetectionMatrix ----

#' Identify Markers and Signatures Present in Individual Samples
#'
#' Create matrices of detection events for individual gene features or combination thereof (i.e., signatures).
#'
#' @rdname makeDetectionMatrices
#'
#' @details
#' The `makeMarkerDetectionMatrix` function returns a marker by sample matrix declaring a feature as detected if it is detected above a given threshold in a specific assay (e.g., `"counts"`, `"logcounts"`).
#'
#' The `makeSignatureDetectionMatrix` function returns a sample by signature matrix declaring a signature (composed of one or more gene features) as detected if all the associated features are detected.
#'
#' The `makeMarkerProportionMatrix` function returns a signature by cluster matrix indicating the proportion of samples expressing detectable levels of each marker in individual predefined clusters.
#'
#' The `makeMarkerProportionScree` function compute the 'cumulative' (i.e., combined) detection rate of markers:
#' the proportion of samples with detectable levels of the first marker, both of the first two markers, etc.
#'
#' @param se An object of class inheriting from [`SummarizedExperiment`][RangedSummarizedExperiment-class].
#' @param markers A character vector, subset of `rownames(se)`.
#' @param threshold Value _above which_ the marker is considered detected.
#' @param assay.type A string specifying which assay values to use, e.g., `"counts"` or `"logcounts"`.
#'
#' @return
#' \describe{
#' \item{`makeMarkerDetectionMatrix`}{A `logical` matrix indicating detectable levels of each marker in each sample.}
#' \item{`makeSignatureDetectionMatrix`}{A `logical` matrix indicating detectable levels of each signature in each sample.}
#' \item{`makeMarkerProportionMatrix`}{A matrix indicating for each feature the proportion of samples expressing detectable levels in each cluster.}
#' \item{`makeMarkerProportionScree`}{A `double` vector indicating the proportion of samples positive for markers in the input `logical` matrix.}
#' }
#'
#' @export
#' @importFrom SummarizedExperiment assay
#'
#' @author Kevin Rue-Albrecht, with C++ code by Aaron Lun
#'
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' nsamples <- 100
#' u <- matrix(rpois(20000, 2), ncol=nsamples)
#' rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
#' se <- SummarizedExperiment(assays=list(counts=u))
#'
#' geneLists <- list(
#'   "Cell type 1" = c("Gene001", "Gene002"),
#'   "Cell type 2" = c("Gene003", "Gene004")
#' )
#' bs <- as(geneLists, "Sets")
#'
#' # Example usage ----
#'
#' markerMatrix <- makeMarkerDetectionMatrix(se, ids(elementInfo(bs)))
#' signatureMatrix <- makeSignatureDetectionMatrix(markerMatrix, bs)
#'
#' tab <- makeMarkerProportionScree(markerMatrix)
#' plot(tab$markers, tab$proportion, ylim=c(0, 1), main="Combined detection scree")
#'
#' se$cluster <- factor(sample(c("A", "B"), ncol(se), TRUE))
#' proportionMatrix <- makeMarkerProportionMatrix(
#'     se, cluster.col="cluster", assay.type="counts", threshold=0)
#' head(proportionMatrix)
makeMarkerDetectionMatrix <- function(
    se, markers, threshold=0, assay.type="counts"
) {
    if (any(duplicated(markers))) {
        warning("Dropping duplicated markers values")
        markers <- unique(markers)
    }
    # Subset the requested assay to the markers of interest
    assayMatrix <- assay(se, assay.type)[markers, , drop=FALSE]
    markerDetectionMatrix <- (assayMatrix > threshold)
    markerDetectionMatrix
}

# makeMarkerProportionScree ----

#' @rdname makeDetectionMatrices
#'
#' @export
makeMarkerProportionScree <- function(matrix) {
    topDetected <- .Call(cxx_num_detected_markers, matrix, seq_len(nrow(matrix))-1L, 1L)
    tab <- rle(sort(topDetected, decreasing=TRUE))
    outTable <- data.frame(
        proportion=cumsum(tab$lengths / ncol(matrix)),
        markers=tab$values)
    outTable
}

# makeSignatureDetectionMatrix ----

#' @rdname makeDetectionMatrices
#'
#' @param matrix A logical matrix indicating the presence of each marker (row) in each sample (column).
#' @param object A collection of signatures inheriting from "[`GeneSetCollection`]" or "[`Sets`]".
#'
#' @export
#'
#' @importFrom S4Vectors FilterRules evalSeparately
makeSignatureDetectionMatrix <- function(
    matrix, object
) {
    filterExpressions <- makeFilterExpression(object)
    fr <- FilterRules(filterExpressions)
    es <- evalSeparately(fr, as.data.frame(t(as.matrix(matrix))))
    es
}

# makeMarkerProportionMatrix ----

#' @rdname makeDetectionMatrices
#'
#' @param cluster.col Name of a column in `colData(se)` that contains
#' a factor indicating cluster membership for each column (i.e. sample) in `se`.
#'
#' @export
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom BiocGenerics rowSums
#' @importFrom Matrix rowSums
makeMarkerProportionMatrix <- function(
    se, cluster.col, assay.type="counts", threshold=0
) {
    stopifnot(!missing(cluster.col))
    stopifnot(is.factor(colData(se)[, cluster.col, drop=TRUE])) # keep as-is to raise an error referring to variables known to the user
    clusterData <- colData(se)[, cluster.col, drop=TRUE]

    # Compute the proportion of samples in each cluster expressing each marker
    markerDetectionMatrix <- makeMarkerDetectionMatrix(se, rownames(se), threshold, assay.type)

    clusterNames <- levels(clusterData)
    numberCellsInCluster <- table(clusterData)

    proportionPositiveByCluster <- matrix(
        data=NA_real_,
        nrow=nrow(markerDetectionMatrix),
        ncol=length(clusterNames),
        dimnames=list(feature=rownames(markerDetectionMatrix), cluster=clusterNames))
    assayMatrix <- assay(se, assay.type)
    for (clusterName in clusterNames) {
        clusterSamples <- which(colData(se)[, cluster.col] == clusterName)
        detectionMatrix <- assayMatrix[, clusterSamples] > threshold
        nDetected <- rowSums(detectionMatrix)
        proportionPositiveByCluster[, clusterName] <- nDetected / length(clusterSamples)
    }

    proportionPositiveByCluster
}
