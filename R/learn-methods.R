
# learnSignatures ----

#' Method to Learn Signatures from SummarizedExperiment
#'
#' These method signatures learn gene set signatures optionally augmented with
#' (semi-)quantitative information for the prediction of sample and cell identities
#' in \code{SummarizedExperiment} objects.
#'
#' @rdname learnHancock
#'
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#' @param method Learning method. See section "Learning methods".
#' @param ... Additional arguments affecting the learning method.
#'
#' @section Learning methods:
#' \describe{
#' \item{PositiveProportionDifference, PPD}{
#' \emph{Requires prior cluster membership information.}
#' Computes the proportion of samples positive for each feature in each cluster.
#' Identifies for each cluster the top \code{n} features showing the maximal difference
#' between the frequency of detection in the cluster of interest and the maximal frequency of detection in any other cluster.}
#' }
#'
#' @return A \code{\link{tbl_geneset}}.
#'
#' @export
#' @importFrom utils head
#' @importFrom Matrix rowSums
#' @importFrom Biobase rowMax
#'
#' @seealso \code{\link{learnMarkersByPositiveProportionDifference}}
#'
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' nsamples <- 100
#' u <- matrix(rpois(20000, 2), ncol=nsamples)
#' rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
#' colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))
#' se <- SummarizedExperiment(assays=list(counts=u))
#'
#' # Example usage ----
#' se1 <- se
#' colData(se1)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se1), replace=TRUE))
#' learnSignatures(se1, method="PositiveProportionDifference", cluster.col="cluster")
learnSignatures <- function(
    se, assay.type="counts", method=c("PositiveProportionDifference", "PPD"), ...
) {
    method <- match.arg(method)

    if (method %in% c("PositiveProportionDifference", "PPD")) {
        out <- learnMarkersByPositiveProportionDifference(se, assay.type=assay.type, ...)
    }

    out
}

# learnMarkersByPositiveProportionDifference ----

#' Identify Markers by Largest Difference of Detection Rate in Clusters
#'
#' This function computes the detection rate of each feature in each cluster.
#' For each cluster, it ranks all the features by decreasing difference between
#' the detection rate in the target cluster, and the maximal detection rate in any other cluster.
#' The function returns up to \code{n} markers for each cluster.
#'
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param cluster.col Name of a column in \code{colData(se)} that contains
#' a factor indicating cluster membership for each column (i.e. sample) in \code{se}.
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#' @param threshold Value \emph{above which} the marker is considered detected.
#' @param n Maximal number of markers allowed for each signature.
#' @param min.diff Minimal difference in detection rate between the target cluster
#' and the summarized detection rate in any other cluster (in the range 0-1).
#' See argument \code{diff.method} below.
#' @param min.prop Minimal proportion of samples in the target cluster where the combined set of markers is detected.
#' @param diff.method Method to contrast the detection rate in the target cluster to that of all other clusters.
#' See Details section.
#'
#' @details
#' \code{diff.method} affects how the detection rate in all clusters \emph{other than the target one} are summarized before comparison with the detection in the target cluster.
#' It is possible to use the minimal (\code{"min"}), \code{"mean"}, \code{"median"} (minimal), or maximal (\code{"max"}) difference between the detection in the target cluster all that of any other cluster.
#'
#' @return A collection of signatures as a "\code{\link{tbl_geneset}}".
#'
#' @export
#' @importFrom Biobase rowMax rowMin
#' @importFrom matrixStats rowMedians
#'
#' @author Kevin Rue-Albrecht
#'
#' @seealso \code{\link{learnSignatures}}
#'
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' nsamples <- 100
#' u <- matrix(rpois(20000, 1), ncol=nsamples)
#' rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
#' colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))
#' se <- SummarizedExperiment(assays=list(counts=u))
#'
#' colData(se)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se), replace=TRUE))
#'
#' # Example usage ----
#' tgs <- learnMarkersByPositiveProportionDifference(se, cluster.col="cluster")
learnMarkersByPositiveProportionDifference <- function(
    se, cluster.col, assay.type="counts", threshold=0, n=Inf, min.diff=0.1, min.prop=0.1,
    diff.method=c("min", "mean", "median", "max")
) {
    # Sanity checks
    if (missing(cluster.col)) {
        stop("cluster.col is required for method 'PositiveProportionDifference'")
    }
    if (min.diff > 1 | min.diff < 0) {
        stop("min.diff must be a scalar in the range [0,1].")
    }
    diff.method <- match.arg(diff.method)
    diffMethods <- c(
        "min"=Biobase::rowMax,# minimal difference: compare to the maximal value
        "mean"=Matrix::rowMeans,
        "median"=matrixStats::rowMedians,
        "max"=Biobase::rowMin # maximal difference: compare to the minimal value
        )
    diff.method <- diffMethods[[diff.method]]

    proportionPositiveByCluster <- makeMarkerProportionMatrix(se, cluster.col, assay.type, threshold)

    getTopMarkers <- function(clusterName, top=n) {
        # Detection rate in the target cluster, and maximum in any other cluster
        df <- data.frame(
            freqTarget=proportionPositiveByCluster[, clusterName, drop=TRUE],
            freqOtherMax=diff.method(proportionPositiveByCluster[
                , setdiff(colnames(proportionPositiveByCluster), clusterName),
                drop=FALSE]),
            row.names=rownames(se)
        )
        # Difference of detection rate
        df$diffFreq <- df$freqTarget - df$freqOtherMax
        # Exclude markers below the minimal difference threshold
        if (!is.na(min.diff)) {
            df <- df[df$diffFreq >= min.diff, , drop=FALSE]
        }
        # 'Combinatorial' detection rate in the target cluster
        # Do not move above the exclusion on min.diff, to save time
        if (!is.na(min.prop)) {
            seSubset <- se[, colData(se)[, cluster.col] == clusterName]
            markerDetectionMatrix <- makeMarkerDetectionMatrix(seSubset, rownames(df), threshold, assay.type)
            # Order markers by decreasing detection rate
            orderedMarkers <- rownames(markerDetectionMatrix)[order(rowSums(markerDetectionMatrix), decreasing=TRUE)]
            markerDetectionMatrix <- markerDetectionMatrix[orderedMarkers, ]
            proportionScreen <- makeMarkerProportionScree(markerDetectionMatrix)
            df <- df[orderedMarkers, , drop=FALSE]
            df$combinedProp <- proportionScreen
            df <- df[df$combinedProp >= min.prop, , drop=FALSE]
        }
        # Reorder the remaining markers.
        # Do not move higher above, it saves time.
        df <- df[order(df$diffFreq, decreasing=TRUE), , drop=FALSE]
        # Extract the request number of markers (default: all)
        head(rownames(df), top)
    }

    clusterNames <- colnames(proportionPositiveByCluster)
    markers <- lapply(clusterNames, getTopMarkers)
    names(markers) <- clusterNames

    tbl <- do.call(tbl_geneset, markers)
    tbl
}

