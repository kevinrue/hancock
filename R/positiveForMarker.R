
#' @rdname positiveForMarker
#' @importFrom SummarizedExperiment assay
setMethod(
    "positiveForMarker", c("SummarizedExperiment"),
    function(x, row, threshold=0, assay="counts"){
        x <- assay(x, assay)
        
        # Call internal function
        return(.positiveForMarker(x, row, threshold))
    }
)

#' @rdname positiveForMarker
setMethod(
    "positiveForMarker", c("matrix"),
    function(x, row, threshold=0){
        # Call internal function
        return(.positiveForMarker(x, row, threshold))
    }
)

#' Detect cells positive for a given marker
#' 
#' @rdname INTERNAL_positiveForMarker
#'
#' @param matrix An expression matrix.
#' @param row Row index of a marker.
#' @param threshold Value \emph{above which} the marker is considered detected.
#'
#' @return A logical vector where \code{TRUE} indicates detection of the marker.
.positiveForMarker <- function(matrix, row, threshold) {
    if (!is.numeric(row)) { row <- match(row, rownames(matrix)) }
    
    isCellPositive <- (matrix[row, ] > threshold)
    
    return(isCellPositive)
}