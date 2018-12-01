
#' @rdname positiveForMarker
#' @importFrom SummarizedExperiment assay
setMethod(
    "positiveForMarker", c("SummarizedExperiment"),
    function(x, row, threshold=0, assay="counts"){
        # Extract assay matrix
        x <- assay(x, assay)
        isCellPositive <- (x[row, ] > threshold)
        isCellPositive
    }
)

#' Identify Markers and Signatures Present in Individual Samples
#'
#' @rdname makeDetectionMatrices
#'
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param markers A character vector, subset of
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
    markerDetectionMatrix <- matrix(
        data=FALSE,
        nrow=ncol(se), ncol=length(markers),
        dimnames=list(colnames(se), markers))
    # TODO: lapply, BiocParallel, ...
    for (gene in markers) {
        markerDetectionMatrix[, gene] <- positiveForMarker(se, gene, threshold, assay.type)
    }
    markerDetectionMatrix
}

#' @rdname makeDetectionMatrices
#'
#' @param matrix A logical matrix indicating the presence of each marker in each sample.
#' See \code{\link{makeMarkerDetectionMatrix}}
#' @param object A set of signatures of class inheriting from "\code{\link{GeneSetCollection}}".
#'
#' @export
#'
#' @importFrom S4Vectors FilterRules evalSeparately
makeSignatureDetectionMatrix <- function(matrix, object) {
    filterExpressions <- .makeFilterExpressionFromGeneSetCollection(object)
    fr <- FilterRules(filterExpressions)
    es <- evalSeparately(fr, as.data.frame(matrix))
    es
}
