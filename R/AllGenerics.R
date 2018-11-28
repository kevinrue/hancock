# positiveForMarker ----

#' Detect cells positive for a given marker
#'
#' @rdname positiveForMarker
#'
#' @param x A \code{\link{SummarizedExperiment}} or \code{matrix}.
#' @param row Row index of a marker.
#' @param threshold Value \emph{above which} the marker is considered detected.
#' @param ... Additional arguments passed to and from methods.
#' @param assay Name of the \code{assay} slot to use.
#'
#' @return A logical vector where \code{TRUE} indicates detection of the marker.
#' @export
#'
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' ncells <- 100
#' u <- matrix(rpois(20000, 2), ncol=ncells)
#' sce <- SummarizedExperiment(assays=list(counts=u))
#'
#' # Example usage ----
#' positiveForMarker(sce, 1, 0)
setGeneric(
    "positiveForMarker", signature = c("x"),
    function(x, row, threshold=0, ...)
        standardGeneric("positiveForMarker")
)
