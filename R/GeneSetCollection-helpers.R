
# uniqueMarkers ----

#' @rdname uniqueMarkers
#' @importFrom GSEABase geneIds
setMethod(
    "uniqueMarkers", c("GeneSetCollection"),
    function(object){
        # NOTE: later, we may trim gene sets to features present in `se`
        # NOTE: in which case, gene sets trimmed to length 0 would have to be dropped (!)
        uniqueMarkers <- unique(unlist(geneIds(object)))
        uniqueMarkers
    }
)

# uniqueSetNames ----

#' @rdname uniqueSetNames
setMethod(
    "uniqueSetNames", c("GeneSetCollection"),
    function(object){
        uniqueSetNames <- names(object)
        uniqueSetNames
    }
)

# makeFilterExpression ----

#' @rdname makeFilterExpression
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
