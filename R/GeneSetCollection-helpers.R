
#' Extract Unique Markers From a GeneSetCollection
#'
#' @rdname INTERNAL_uniqueMarkers
#'
#' @param object A GeneSetCollection object.
#'
#' @return A character vector of unique markers across all signatures.
#'
#' @author Kevin Rue-Albrecht
.uniqueMarkers <- function(
    object
) {
    # NOTE: later, we may trim signatures to features present in `se`
    # NOTE: in which case, signatures trimmed to length 0 would have to be dropped (!)
    uniqueMarkers <- unique(unlist(geneIds(object)))
    uniqueMarkers
}

#' Build Filter Expressions from GeneSetCollection objects
#'
#' @rdname INTERNAL_makeFilterExpressionFromGeneSetCollection
#'
#' @param object A \code{\link{GeneSetCollection}}.
#'
#' @return A list of \code{\link{expression}} that combines the markers
#' listed in each signature.
#'
#' @author Kevin Rue-Albrecht
.makeFilterExpressionFromGeneSetCollection <- function(object) {
    filterExpressions <- lapply(
        names(object),
        function(x){ parse(text=paste(geneIds(object[[x]]), collapse=" & ")) }
    )
    names(filterExpressions) <- names(object)
    filterExpressions
}
