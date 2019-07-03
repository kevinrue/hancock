
# uniqueMarkerNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueMarkerNames,Sets-methods
#'
#' @export
#' @importFrom unisets elementInfo ids
setMethod(
    "uniqueMarkerNames", c("Sets"),
    function(object){
        # NOTE: later, we may trim gene sets to features present in `se`
        # NOTE: in which case, gene sets trimmed to length 0 would have to be dropped (!)
        uniqueMarkerNames <- ids(elementInfo(object))
        uniqueMarkerNames
    }
)

# uniqueSetNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueSetNames,Sets-methods
#'
#' @export
#' @importFrom unisets setInfo ids
setMethod(
    "uniqueSetNames", c("Sets"),
    function(object){
        uniqueSetNames <- ids(setInfo(object))
        uniqueSetNames
    }
)

# makeFilterExpression ----

#' @rdname makeFilterExpression
#' @aliases makeFilterExpression,Sets-methods
#'
#' @export
#' @importFrom unisets setInfo ids
setMethod(
    "makeFilterExpression", c("Sets"), function(object){

        xList <- as(object, "list")

        .buildSingleExpression <- function(setName) {
            geneIds <- ids(xList[[setName]])
            parse(text=paste(sprintf("`%s`", geneIds), collapse=" & "))
        }

        filterExpressions <- lapply(
            ids(setInfo(object)),
            .buildSingleExpression
        )
        names(filterExpressions) <- ids(setInfo(object))
        filterExpressions
    }
)
