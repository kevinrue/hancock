
# uniqueMarkerNames ----

#' Extract Unique Signature Names from a Collection of Gene Sets
#'
#' `uniqueMarkerNames` and `uniqueSetNames` extract the character vectors of unique marker and set names from objects that store collections of gene sets.
#'
#' @rdname uniqueMarkerNames
#' @aliases uniqueMarkerNames
#'
#' @param object An object of class inheriting from [`BaseSets-class`] or [`GeneSetCollection-class`].
#'
#' @return A character vector of unique set or marker names across all gene sets.
#' @export
#'
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#'
#' library(unisets)
#' gmt <- system.file(package = "hancock", "extdata", "seurat_pbmc3k.gmt")
#' genesets <- import(gmt)
#'
#' # Usage ----
#' um <- uniqueMarkerNames(genesets)
#' us <- uniqueSetNames(genesets)
setGeneric(
    "uniqueMarkerNames", signature=c("object"),
    function(object)
        standardGeneric("uniqueMarkerNames")
)

# uniqueSetNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueSetNames
#'
#' @export
setGeneric(
    "uniqueSetNames", signature=c("object"),
    function(object)
        standardGeneric("uniqueSetNames")
)

# makeFilterExpression ----

#' Build Filter Expressions from `GeneSetCollection` objects
#'
#' This function create a list of unevaluated expressions representing a collection of signatures.
#' The resulting expressions can be evaluated as `FilterRules` inside an environment such as a `data.frame` of sample-by-gene detection events.
#'
#' @rdname makeFilterExpression
#'
#' @param object An object of class inheriting from [`BaseSets`] or [`GeneSetCollection`].
#'
#' @return A list of [`expression`] that combines the markers listed in each gene set.
#'
#' @export
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#'
#' library(unisets)
#' gmt <- system.file(package = "hancock", "extdata", "seurat_pbmc3k.gmt")
#' genesets <- import(gmt)
#'
#' # Example usage ----
#' fe <- makeFilterExpression(genesets)
setGeneric(
    "makeFilterExpression", signature=c("object"),
    function(object)
        standardGeneric("makeFilterExpression")
)
