
# uniqueMarkerNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueMarkerNames,GeneSetCollection-methods
#'
#' @export
#' @importFrom GSEABase geneIds
setMethod(
    "uniqueMarkerNames", c("GeneSetCollection"),
    function(object){
        # NOTE: later, we may trim gene sets to features present in `se`
        # NOTE: in which case, gene sets trimmed to length 0 would have to be dropped (!)
        uniqueMarkerNames <- unique(unlist(geneIds(object)))
        uniqueMarkerNames
    }
)

# uniqueSetNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueSetNames,GeneSetCollection-methods
#'
#' @export
setMethod(
    "uniqueSetNames", c("GeneSetCollection"),
    function(object){
        uniqueSetNames <- names(object)
        uniqueSetNames
    }
)

# makeFilterExpression ----

#' @rdname makeFilterExpression
#' @aliases makeFilterExpression,GeneSetCollection-methods
#'
#' @export
#' @importFrom GSEABase geneIds
setMethod(
    "makeFilterExpression", c("GeneSetCollection"),
    function(object){

        .buildSingleExpression <- function(setName) {
            parse(text=paste(sprintf("`%s`", geneIds(object[[setName]])), collapse=" & "))
        }

        filterExpressions <- lapply(
            names(object),
            .buildSingleExpression
        )
        names(filterExpressions) <- names(object)
        filterExpressions
    }
)
