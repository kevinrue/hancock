
#' @describeIn predictProportionSignatureByCluster
#' Visualize proportion of samples positive signatures in each cluster.
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
