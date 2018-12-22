
# learnSignaturess ----

#' Method to Learn Signatures from SummarizedExperiment
#'
#' These method signatures learn gene set signatures optionally augmented with
#' (semi-)quantitative information for the prediction of sample and cell identities
#' in \code{SummarizedExperiment} objects.
#'
#' @rdname learnSignatures
#'
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#' @param method Learning method. See section "Learning methods".
#' @param ... Additional arguments affecting the learning method.
#'
#' @section Learning methods:
#' \describe{
#' \item{ProportionDifference}{
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
#' colData(se1)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se1), replace = TRUE))
#' learnSignatures(se1, method="ProportionDifference", cluster.col="cluster")
learnSignatures <- function(
    se, assay.type="counts", method=c("ProportionDifference"), ...
) {
    method <- match.arg(method)

    if (identical(method, "ProportionDifference")) {
        out <- learnSignaturesByProportionDifference(se, assay.type=assay.type, ...)
    }

    out
}

# learnSignaturesByProportionDifference ----

#' Identify Markers by Largest Difference of Detection Rate in Clusters
#'
#' This function computes the detection rate of each gene feature in each cluster.
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
#'
#' @return A collection of signatures as a "\code{\link{tbl_geneset}}".
#'
#' @export
#'
#' @author Kevin Rue-Albrecht
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' nsamples <- 100
#' u <- matrix(rpois(20000, 1), ncol=nsamples)
#' rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
#' colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))
#' se <- SummarizedExperiment(assays=list(counts=u))
#'
#' colData(se)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se), replace = TRUE))
#'
#' # Example usage ----
#' tgs <- learnSignaturesByProportionDifference(se, cluster.col="cluster")
learnSignaturesByProportionDifference <- function(
    se, cluster.col, assay.type="counts", threshold=0, n=2
) {
    # Sanity checks
    if (missing(cluster.col)) {
        stop("cluster.col is required for method 'ProportionDifference'")
    }
    stopifnot(!missing(cluster.col))
    stopifnot(is.factor(colData(se)[, cluster.col, drop=TRUE]))
    clusterData <- colData(se)[, cluster.col, drop=TRUE]

    # Compute the proportion of each cluster positive for each marker
    markerDetectionMatrix <- makeMarkerDetectionMatrix(se, rownames(se), threshold, assay.type)

    clusterNames <- levels(clusterData)
    numberCellsInCluster <- table(clusterData)

    proportionPositiveByCluster <- matrix(
        data=NA_real_,
        nrow=nrow(markerDetectionMatrix),
        ncol=length(clusterNames),
        dimnames=list(feature=rownames(markerDetectionMatrix), cluster=clusterNames))
    x <- assay(se, assay.type)
    for (clusterName in clusterNames) {
        clusterSamples <- which(colData(se)[, cluster.col] == clusterName)
        nDetected <- Matrix::rowSums(x[, clusterSamples] > threshold)
        proportionPositiveByCluster[, clusterName] <- nDetected / length(clusterSamples)
    }

    markers <- lapply(clusterNames, function(clusterName, top=n) {
        freqTarget <- proportionPositiveByCluster[, clusterName]
        freqMaxOther <- rowMax(proportionPositiveByCluster[, setdiff(colnames(proportionPositiveByCluster), clusterName)])
        diffFreq <- freqTarget - freqMaxOther
        rankedGenes <- rownames(se)[order(diffFreq, decreasing = TRUE)]
        head(rankedGenes, top)
    })
    names(markers) <- clusterNames

    tbl <- do.call(tbl_geneset, markers)
    tbl
}

