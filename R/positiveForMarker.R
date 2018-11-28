
#' @rdname positiveForMarker
#' @importFrom SummarizedExperiment assay
setMethod(
    "positiveForMarker", c("SummarizedExperiment"),
    function(x, row, threshold=0, assay="counts"){
        # Extract assay matrix
        x <- assay(x, assay)
        isCellPositive <- (x[row, ] > threshold)
        isCellPositive
    }
)
