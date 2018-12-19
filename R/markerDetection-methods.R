
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
#' ncells <- 100
#' u <- matrix(rpois(20000, 2), ncol=ncells)
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
makeSignatureDetectionMatrix <- function(matrix, object) {
    filterExpressions <- makeFilterExpression(object)
    fr <- FilterRules(filterExpressions)
    es <- evalSeparately(fr, as.data.frame(t(matrix)))
    es
}
