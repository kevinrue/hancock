
# predict ----

#' Predict Method for Collection of Gene Set Signatures
#'
#' These method signatures apply gene set signatures optionally augmented with
#' (semi-)quantitative information to the prediction of sample and cell identities
#' in [`SummarizedExperiment`][RangedSummarizedExperiment-class] objects.
#'
#' @name predict-hancock
#' @rdname predictSignatures
#' @aliases predict predict-methods
#'
#' @param object A set of signatures of class inheriting from [`BaseSets-class`] or [`GeneSetCollection-class`].
#' @param se An object of class inheriting from [`SummarizedExperiment`][RangedSummarizedExperiment-class].
#' @param assay.type A string specifying which assay values to use, e.g., `"counts"` or `"logcounts"`.
#' @param method Prediction method. See section "Prediction methods".
#' @param ... Additional arguments affecting the predictions produced.
#'
#' @section Prediction methods:
#' \describe{
#' \item{ProportionPositive, PP}{
#' _Requires prior cluster membership information._
#' Computes the proportion of samples positive for each signature in each cluster.
#' Assigns to each cluster the signature detected in the highest proportion of samples.}
#' }
#'
#' @return The object `se`, updated as follows:
#' \itemize{
#' \item in the `metadata` slot, a `"hancock"` item is added (or updated) with information tracing the prediction process (e.g., method, signatures).
#' \item in the `"colData"` slot, a `DataFrame` is nested in a new (or updated) `"hancock"` column.
#' This DataFrame contains predicted labels in the first column and additional information in further columns for each column in `se`.
#' }
#'
#' @export
#' @importFrom S4Vectors metadata
#'
#' @seealso [`predictByProportionPositive`].
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
#' bs <- BaseSets(
#'     relations=DataFrame(
#'         element = c("Gene001", "Gene002", "Gene003", "Gene004"),
#'         set     = c(rep("Cell type 1", 2), rep("Cell type 2", 2))
#'     )
#' )
#'
#' # Example usage ----
#'
#' se1 <- se
#' colData(se1)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se1), replace=TRUE))
#' se1 <- predict(bs, se1, method="ProportionPositive", cluster.col="cluster")
#' # Visualise the count of samples predicted for each signature in each cluster
#' barplotPredictionCount(se1, highlight=c("Cell type 1"))
#'
#' barplotPredictionProportion(se1, highlight=c("Cell type 2"))
#'
#' library(SingleCellExperiment)
#' sce1 <- as(se1, "SingleCellExperiment")
#' reducedDim(sce1, "PCA") <- prcomp(t(assay(sce1)))$x
#' reducedDimPrediction(sce1, highlight="Cell type 1", redDimType="PCA", x=1, y=2)
predict.GeneSetCollection <- function(
    object, se, assay.type="counts", method=c("ProportionPositive", "PP"), ...
) {
    .predictAnyGeneSetClass(object, se, assay.type, method, ...)
}

#' @rdname predictSignatures
#'
#' @export
#'
#' @method predict tbl_geneset
predict.tbl_geneset <- function(
    object, se, assay.type="counts", method=c("ProportionPositive", "PP"), ...
) {
    .predictAnyGeneSetClass(object, se, assay.type, method, ...)
}

#' @rdname predictSignatures
#' @export
predict.BaseSets <- function(
    object, se, assay.type="counts", method=c("ProportionPositive", "PP"), ...
) {
    .predictAnyGeneSetClass(object, se, assay.type, method, ...)
}

#' Internal Predict Method for Any Type of Gene Set Collection
#'
#' This function is called by both [`BaseSets-class`] or [`GeneSetCollection-class`] signatures.
#' Dispatch occurs in downstream functions (e.g., [`uniqueMarkerNames()`], [`uniqueSetNames()`]).
#'
#' @param object A collection of signatures inheriting from [`BaseSets-class`] or [`GeneSetCollection-class`].
#' @param se An object of class inheriting from [`SummarizedExperiment`][RangedSummarizedExperiment-class].
#' @param assay.type A string specifying which assay values to use, e.g., `"counts"` or `"logcounts"`.
#' @param method Prediction method. See section "Prediction methods".
#' @param ... Additional arguments affecting the predictions produced.
#'
#' @rdname INTERNAL_predictAnyGeneSetClass
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom utils packageVersion
#'
#' @author Kevin Rue-Albrecht
.predictAnyGeneSetClass <- function(
    object, se, assay.type="counts", method=c("ProportionPositive", "PP"), ...
){
    method <- match.arg(method)

    # NOTE: match.arg above ensures that invalid methods throw an error
    if (method %in% c("ProportionPositive", "PP")) {
        se <- predictByProportionPositive(object, se, ..., assay.type=assay.type)
    }

    # Update the hancock metadata.
    # Add global new (general) metadata before existing (method-specific) metadata
    existingMetadata <- metadata(se)[[getPackageName()]]
    newMetadata <- list(
        GeneSets=object,
        method=method,
        packageVersion=packageVersion(getPackageName())
    )
    metadata(se)[[getPackageName()]] <- append(newMetadata, existingMetadata)

    se
}

# predictByProportionPositive ----

#' Identify the Dominant Signatures in Clusters of Samples
#'
#' The [`predictByProportionPositive()`] function computes the proportion of samples positive for each signature in each (predefined) cluster
#' and identifies the predominant signature in each cluster.
#' The function stores information tracing the prediction process in the `metadata` slot. See Details.
#'
#' @name predictByProportionPositive
#' @rdname predictByProportionPositive
#'
#' @details
#' The function populates the `"hancock"` element of the `metadata` slot with the following fields and values:
#' \describe{
#' \item{`"GeneSets"`}{Signatures used to make the predictions}
#' \item{`"method"`}{Name of the method used to make the predictions}
#' \item{`"packageVersion"`}{`hancock` version used to make the predictions}
#' \item{`"ProportionPositiveByCluster"`}{Matrix indicating the proportion of samples in each cluster that are positive for each signature.}
#' \item{`"TopSignatureByCluster"`}{Named vector indicating the predominant signature for each cluster.}
#' }
#'
#' @param object A collection of signatures inheriting from [`BaseSets-class`] or [`GeneSetCollection-class`]".
#' @param se An object of class inheriting from [`SummarizedExperiment`][RangedSummarizedExperiment-class].
#' @param cluster.col Name of a column in `colData(se)` that contains
#' a factor indicating cluster membership for each column (i.e. sample) in `se`.
#' @param assay.type A string specifying which assay values to use, e.g., `"counts"` or `"logcounts"`.
#' @param threshold Value _above which_ the marker is considered detected.
#'
#' @return The object `se`, updated as follows:
#' \itemize{
#' \item in the `metadata` slot, a `"hancock"` item is added (or updated) with information tracing the prediction process. See Details.
#' \item in the `colData` slot, a `DataFrame` is nested in a new (or updated) `"hancock"` column.
#' This DataFrame contains predicted labels in the first and only column.
#' }
#'
#' @importFrom SummarizedExperiment colData colData<- metadata<-
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
#' @author Kevin Rue-Albrecht
#'
#' @seealso [`predict.GeneSetCollection()`], [`predict.BaseSets()`].
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
#' bs <- BaseSets(
#'     relations=DataFrame(
#'         element = c("Gene001", "Gene002", "Gene003", "Gene004"),
#'         set     = c(rep("Cell type 1", 2), rep("Cell type 2", 2))
#'     )
#' )
#' colData(se)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se), replace=TRUE))
#'
#' # Example usage ----
#' library(circlize)
#' # Identify the dominant signature in each cluster
#' se <- predictByProportionPositive(bs, se, cluster.col="cluster")
#' # Visualise the proportion of samples positive for each signature in each cluster
#' plotProportionPositive(
#'   se, cluster_rows=FALSE, cluster_columns=FALSE,
#'   col=colorRamp2(c(0, 100), c("white", "red")))
predictByProportionPositive <- function(
    object, se, cluster.col, assay.type="counts", threshold=0
) {
    # Sanity checks
    if (missing(cluster.col)) {
        stop("cluster.col is required for method 'ProportionPositive'")
    }
    stopifnot(!missing(cluster.col))
    stopifnot(is.factor(colData(se)[, cluster.col, drop=TRUE]))
    clusterData <- colData(se)[, cluster.col, drop=TRUE]

    # Compute the proportion of each cluster positive for each signature
    uniqueMarkersIds <- uniqueMarkerNames(object)
    stopifnot(all(uniqueMarkersIds %in% rownames(se)))
    markerDetectionMatrix <- makeMarkerDetectionMatrix(se, uniqueMarkersIds, threshold, assay.type)
    signatureMatrix <- makeSignatureDetectionMatrix(markerDetectionMatrix, object)

    clusterNames <- levels(clusterData)
    numberCellsInCluster <- table(clusterData)
    signatureNames <- uniqueSetNames(object)

    proportionPositiveByCluster <- matrix(
        data=0,
        nrow=length(clusterNames),
        ncol=ncol(signatureMatrix),
        dimnames=list(cluster=clusterNames, signature=colnames(signatureMatrix)))
    for (signatureName in colnames(signatureMatrix)) {
        countSignatureTable <- table(clusterData, signatureMatrix[, signatureName])
        if ("TRUE" %in% colnames(countSignatureTable)) {
            proportionPositiveByCluster[, signatureName] <- countSignatureTable[, "TRUE"] / numberCellsInCluster
        }
    }

    # For each cluster, identify the most frequent signature
    # TODO: warning if ties
    maxColIdx <- max.col(proportionPositiveByCluster, ties.method="first")
    maxSignatureName <- factor(
        x=colnames(proportionPositiveByCluster)[maxColIdx],
        levels=signatureNames)
    names(maxSignatureName) <- rownames(proportionPositiveByCluster)

    # Assign most frequent signature to every cell in each cluster
    newColData <- DataFrame(
        prediction=maxSignatureName[colData(se)[, cluster.col, drop=TRUE]]
    )
    colData(se)[[getPackageName()]] <- newColData

    # Store the proportion of cluster positive for each signature in metadata for later plotting
    newMetadata <- list(
        ProportionPositiveByCluster=proportionPositiveByCluster,
        TopSignatureByCluster=maxSignatureName
    )
    metadata(se)[[getPackageName()]] <- newMetadata

    se
}
