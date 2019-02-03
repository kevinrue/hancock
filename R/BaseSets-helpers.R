
# uniqueMarkerNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueMarkerNames
#'
#' @export
#' @importFrom unisets elementData ids
setMethod(
    "uniqueMarkerNames", c("BaseSets"),
    function(object){
        # NOTE: later, we may trim gene sets to features present in `se`
        # NOTE: in which case, gene sets trimmed to length 0 would have to be dropped (!)
        uniqueMarkerNames <- ids(elementData(object))
        uniqueMarkerNames
    }
)

# uniqueSetNames ----

#' @rdname uniqueMarkerNames
#' @aliases uniqueSetNames
#'
#' @export
#' @importFrom unisets setData ids
setMethod(
    "uniqueSetNames", c("BaseSets"),
    function(object){
        uniqueSetNames <- ids(setData(object))
        uniqueSetNames
    }
)

# makeFilterExpression ----

#' @rdname makeFilterExpression
#' @importFrom unisets setData ids
setMethod(
    "makeFilterExpression", c("BaseSets"), function(object){

        xList <- as(object, "list")

        .buildSingleExpression <- function(setName) {
            geneIds <- ids(xList[[setName]])
            parse(text=paste(sprintf("`%s`", geneIds), collapse=" & "))
        }

        filterExpressions <- lapply(
            setIds(object),
            .buildSingleExpression
        )
        names(filterExpressions) <- setIds(object)
        filterExpressions
    }
)
