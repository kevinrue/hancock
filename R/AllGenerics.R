# positiveForMarker ----

#' Detect Samples Positive for Given Markers
#'
#' Declare whether individual gene features (i.e., "markers") are detected above a given threshold in a given assay in each sample.
#'
#' @rdname positiveForMarker
#'
#' @param x A \code{\link{SummarizedExperiment}}.
#' @param row Row index of markers to examine.
#' @param threshold Value \emph{above which} markers are considered detected.
#' @param ... Additional arguments passed to and from methods.
#' @param assay Name of the \code{assay} slot to use.
#'
#' @return A logical vector (single marker) or matrix (multiple markers) where \code{TRUE} indicates detection of the marker.
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
    "positiveForMarker", signature=c("x"),
    function(x, row, threshold=0, ...)
        standardGeneric("positiveForMarker")
)

# uniqueMarkers ----

#' Extract Unique Markers from a Collection of Gene Sets
#'
#' @rdname uniqueMarkers
#'
#' @param object A \code{\link{tbl_geneset}} or \code{GeneSetCollection}.
#'
#' @return A character vector of unique markers across all gene sets.
#' @export
#'
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#' library(GeneSet)
#' tgs <- tbl_geneset(
#'     "Cell type 1" = c("Gene001", "Gene002"),
#'     "Cell type 2" = c("Gene002", "Gene003", "Gene004")
#' )
#'
#' # Example usage ----
#' uniqueMarkers(tgs)
setGeneric(
    "uniqueMarkers", signature=c("object"),
    function(object)
        standardGeneric("uniqueMarkers")
)

# makeFilterExpression ----

#' Build Filter Expressions from GeneSetCollection objects
#'
#' This function create a list of unevaluated expressions representing a collection of signatures.
#' The resulting expressions can be evaluated as \code{FilterRules} inside an environment such as a \code{data.frame} of sample-by-gene detection events.
#'
#' @rdname makeFilterExpression
#'
#' @param object A \code{\link{tbl_geneset}} or \code{\link{GeneSetCollection}}.
#'
#' @return A list of \code{\link{expression}} that combines the markers
#' listed in each gene set.
#'
#' @export
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#' library(GeneSet)
#' tgs <- tbl_geneset(
#'     "Cell type 1" = c("Gene001", "Gene002"),
#'     "Cell type 2" = c("Gene002", "Gene003", "Gene004")
#' )
#'
#' # Example usage ----
#' makeFilterExpression(tgs)
setGeneric(
    "makeFilterExpression", signature=c("object"),
    function(object)
        standardGeneric("makeFilterExpression")
)
