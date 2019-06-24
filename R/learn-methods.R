
# learnSignatures ----

#' Method to Learn Signatures from SummarizedExperiment
#'
#' These method signatures learn gene set signatures optionally augmented with
#' (semi-)quantitative information for the prediction of sample and cell identities
#' in `SummarizedExperiment` objects.
#'
#' @rdname learnSignatures
#'
#' @param se An object of class inheriting from "[`SummarizedExperiment`]".
#' @param assay.type A string specifying which assay values to use, e.g., `"counts"` or `"logcounts"`.
#' @param method Learning method. See section "Learning methods".
#' @param ... Additional arguments affecting the learning method.
#'
#' @section Learning methods:
#' \describe{
#' \item{PositiveProportionDifference, PPD}{
#' _Requires prior cluster membership information._
#' This method computes the proportion of samples positive for each feature in each cluster,
#' and subsequently identifies for each cluster the features showing the maximal difference
#' between the detection rate in the cluster of interest and the detection rate in all other clusters.}
#' }
#'
#' @return A [`Sets`] object.
#'
#' @export
#'
#' @seealso [`learnMarkersByPositiveProportionDifference()`]
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
#' gs <- learnSignatures(se1, method="PositiveProportionDifference", cluster.col="cluster")
#'
#' relations(gs)
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
#' the detection rate in the target cluster, and the detection rate in all other clusters.
#' The function can limit results up to `n` markers for each cluster.
#'
#' @param se An object of class inheriting from "[`SummarizedExperiment`]".
#' @param cluster.col Name of a column in `colData(se)` that contains
#' a factor indicating cluster membership for each column (i.e. sample) in `se`.
#' @param assay.type A string specifying which assay values to use, e.g., `"counts"` or `"logcounts"`.
#' @param threshold Value _above which_ the marker is considered detected.
#' @param n Maximal number of markers allowed for each signature.
#' @param min.diff Minimal difference in detection rate between the target cluster
#' and the summarized detection rate in any other cluster (in the range 0-1).
#' See argument `diff.method`.
#' @param min.prop Minimal proportion of samples in the target cluster where the combined set of markers is detected.
#' @param diff.method Method to contrast the detection rate in the target cluster to that of all other clusters.
#' See Details section.
#'
#' @details
#' `diff.method` controls how the detection rate in all clusters _other than the target one_ are summarized before comparison with the detection in the target cluster.
#' It is possible to rank features using the minimal (`"min"`), `"mean"`, `"median"`, or maximal (`"max"`) difference between the detection rate in the target cluster and those of all other clusters.
#'
#' @return A collection of signatures as a "[`Sets`]" object.
#'
#' @export
#' @importFrom Biobase rowMax rowMin
#' @importFrom matrixStats rowMedians
#' @importFrom utils head
#'
#' @author Kevin Rue-Albrecht
#'
#' @seealso [`learnSignatures`].
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
#'
#' baseset <- learnMarkersByPositiveProportionDifference(se, cluster.col="cluster")
#'
#' relations(baseset)
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
    diffMethodsMap <- c(
        "min"=Biobase::rowMax,# minimal difference: compare to the maximal other value
        "mean"=Matrix::rowMeans,
        "median"=matrixStats::rowMedians,
        "max"=Biobase::rowMin # maximal difference: compare to the minimal other value
        )
    diff.FUN <- diffMethodsMap[[diff.method]]
    diffProportionFieldName <- paste0(diff.method, "DifferenceProportion")
    # Compute the proportion of each cluster positive for each marker
    proportionPositiveByCluster <- makeMarkerProportionMatrix(se, cluster.col, assay.type, threshold)
    # Helper function. Consider refactoring.
    computeClusterMarkers <- function(clusterName, top=n) {
        # Detection rate in the target cluster, and maximum in any other cluster
        df <- data.frame(
            ProportionPositive=proportionPositiveByCluster[, clusterName, drop=TRUE],
            ProportionOthers=diff.FUN(proportionPositiveByCluster[
                , setdiff(colnames(proportionPositiveByCluster), clusterName),
                drop=FALSE]),
            row.names=rownames(se)
        )
        # Difference of detection rate
        df[[diffProportionFieldName]] <- df[["ProportionPositive"]] - df[["ProportionOthers"]]
        # Exclude markers below the minimal difference threshold
        if (!is.na(min.diff)) {
            df <- df[df[[diffProportionFieldName]] >= min.diff, , drop=FALSE]
        }
        # 'Combinatorial' detection rate in the target cluster
        # Do not move above the exclusion on min.diff, to save time
        combinedProportion <- NA_real_
        if (min.prop > 0) {
            # Exclude markers detected alone below the threshold
            df <- df[df[["ProportionPositive"]] >= min.prop, , drop=FALSE]
            # Order markers by decreasing detection rate
            df <- df[order(df[["ProportionPositive"]], decreasing=TRUE), ]
            # Compute which of the remaining markers are detected in each sample
            seSubset <- se[, colData(se)[[cluster.col]] == clusterName]
            markerDetectionMatrix <- makeMarkerDetectionMatrix(seSubset, rownames(df), threshold, assay.type)
            # Identify the maximal number of markers simultaneously detected above the threshold
            proportionScree <- makeMarkerProportionScree(markerDetectionMatrix)
            maxRow <- head(which(proportionScree$proportion >= min.prop), 1)
            maxMarkers <- proportionScree[maxRow, "markers", drop=TRUE]
            maxProportion <- proportionScree[maxRow, "proportion", drop=TRUE]
            df <- head(df, maxMarkers)
            combinedProportion <- maxProportion
        }
        # Reorder the remaining markers by differential detection rate
        # Do not move higher above, it saves time.
        df <- df[order(df[[diffProportionFieldName]], decreasing=TRUE), , drop=FALSE]
        # Extract the request number of markers (default: all)
        out <- list(
            table=head(df, top),
            combinedProportion=combinedProportion
        )
        out
    }
    # Compute the markers for each cluster
    clusterNames <- colnames(proportionPositiveByCluster)
    clusterResults <- lapply(clusterNames, "computeClusterMarkers")
    # Extract tables of marker metadata
    markersTables <- lapply(clusterResults, function(x){x[["table"]]})
    # Extract gene sets
    markersList <- lapply(markersTables, "rownames")
    names(markersList) <- clusterNames
    # Extract combine detection rate of each signature
    setDetectionRates <- lapply(clusterResults, function(x){x[["combinedProportion"]]})
    names(setDetectionRates) <- clusterNames
    # Combine tables of marker metadata (relations metadata)
    metadata <- do.call(rbind, markersTables)
    markersTable <- DataFrame(
        element=unlist(markersList, use.names=FALSE),
        set=rep(names(markersList), lengths(markersList)),
        metadata[, c("ProportionPositive", diffProportionFieldName), drop=FALSE]
    )
    # Collect set metadata
    setInfo <- IdVector(unique(names(markersList)))
    mcols(setInfo) <- DataFrame(
        ProportionPositive=unlist(setDetectionRates)[ids(setInfo)]
    )
    # Collect element metadata
    elementInfo <- IdVector(unique(unlist(markersList)))
    detectionMatrix <- as(assay(se[ids(elementInfo), ], assay.type) > threshold, "Matrix")
    mcols(elementInfo) <- DataFrame(
        ProportionPositive=rowSums(detectionMatrix) / ncol(se)
    )
    # Make a Sets
    Sets(relations=markersTable, elementInfo=elementInfo, setInfo=setInfo)
}
