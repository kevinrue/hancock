
#' @describeIn predictProportionSignatureByCluster
#' Returns a \code{Heatmap} displaying the proportion (on a scale from 0 to 100) of samples that are positive for each individual signature in each cluster.
#'
#' @aliases plotProportionPositive
#'
#' @param ... Additional arguments to be passed to methods.
#'
#' @export
#' @importFrom ComplexHeatmap Heatmap
plotProportionPositive <- function(
    se, ...
){
    ppbc <- metadata(se)[["Hancock"]][["ProportionPositiveByCluster"]]
    if (is.null(ppbc)) {
        stop("Method 'ProportionPositive' was not run yet.")
    }
    Heatmap(matrix=t(ppbc*100), name="Proportion (%)", ...)
}
