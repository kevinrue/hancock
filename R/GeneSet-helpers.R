
# uniqueMarkerNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueMarkerNames,tbl_geneset-methods
#'
#' @export
setMethod(
    "uniqueMarkerNames", "tbl_geneset", function(object) {
      # NOTE: later, we may trim gene sets to features present in `se`
      # NOTE: in which case, gene sets trimmed to length 0 would have to be dropped (!)
      uniqueMarkerNames <- unique(object$gene)
      uniqueMarkerNames
    }
)

# uniqueSetNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueSetNames,tbl_geneset-methods
#'
#' @export
setMethod(
    "uniqueSetNames", "tbl_geneset", function(object) {
        uniqueSetNames <- levels(object$set)
        uniqueSetNames
    }
)

# makeFilterExpression ----

#' @rdname makeFilterExpression
#' @aliases makeFilterExpression,tbl_geneset-methods
#'
#' @export
setMethod(
    "makeFilterExpression", c("tbl_geneset"),
    function(object){
        .buildSingleExpression <- function(setName) {
            geneIds <- object[object$set == setName, "gene", drop=TRUE]
            parse(text=paste(sprintf("`%s`", geneIds), collapse=" & "))
        }

        filterExpressions <- lapply(
            levels(object$set),
            .buildSingleExpression
        )
        names(filterExpressions) <- levels(object$set)
        filterExpressions
    }
)
