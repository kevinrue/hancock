
#' Build Filter Expressions from GeneSetCollection objects
#'
#' @rdname INTERNAL_makeFilterExpressionFromGeneSetCollection
#'
#' @param object A \code{\link{GeneSetCollection}}.
#'
#' @return A list of \code{\link{expression}} that combines the markers
#' listed in each signature.
.makeFilterExpressionFromGeneSetCollection <- function(object) {
    filterExpressions <- lapply(
        names(object),
        function(x){ parse(text=paste(geneIds(object[[x]]), collapse=" & ")) }
    )
    names(filterExpressions) <- names(object)
    filterExpressions
}
