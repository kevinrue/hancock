
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
