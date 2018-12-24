
# predict ----

#' Predict Method for Collection of Gene Set Signatures
#'
#' These method signatures apply gene set signatures optionally augmented with
#' (semi-)quantitative information to the prediction of sample and cell identities
#' in \code{SummarizedExperiment} objects.
#'
#' @rdname predictHancock
#' @aliases predict
#'
#' @param object A set of signatures of class inheriting from "\code{\link{GeneSetCollection}}" or "\code{\link{tbl_geneset}}".
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#' @param method Prediction method. See section "Prediction methods".
#' @param ... Additional arguments affecting the predictions produced.
#'
#' @section Prediction methods:
#' \describe{
#' \item{ProportionPositive, PP}{
#' \emph{Requires prior cluster membership information.}
#' Computes the proportion of samples positive for each signature in each cluster.
#' Assigns to each cluster the signature detected in the highest proportion of samples.}
#' }
#'
#' @return The object \code{se}, updated as follows:
#' \itemize{
#' \item in the \code{metadata} slot, a \code{"Hancock"} item is added (or updated) with information tracing the prediction process (e.g., method, signatures).
#' \item in the \code{colData} slot, a \code{DataFrame} is nested in a new (or updated) \code{"Hancock"} column.
#' This DataFrame contains predicted labels in the first column and additional information in further columns for each column in \code{se}.
#' }
#'
#' @export
#' @method predict GeneSetCollection
#' @importFrom S4Vectors metadata
#'
#' @seealso \code{\link{predictByProportionPositive}}
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
#' gsc <- GeneSetCollection(list(
#'     GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
#'     GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
#' ))
#'
#' # Example usage ----
#' se1 <- se
#' colData(se1)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se1), replace=TRUE))
#' predict(gsc, se1, method="ProportionPositive", cluster.col="cluster")
predict.GeneSetCollection <- function(
    object, se, assay.type="counts", method=c("ProportionPositive", "PP"), ...
) {
    .predictAnyGeneSetClass(object, se, assay.type, method, ...)
}

#' @rdname predictHancock
#' @export
#' @method predict tbl_geneset
predict.tbl_geneset <- function(
    object, se, assay.type="counts", method=c("ProportionPositive", "PP"), ...
) {
    .predictAnyGeneSetClass(object, se, assay.type, method, ...)
}

#' Internal Predict Method for Any Type of Gene Set Collection
#'
#' This function is called by both \code{GeneSetCollection} and \code{tbl_geneset} signatures.
#' Dispatch occurs in downstream functions (e.g. \code{positiveForMarker}, \code{makeSignatureDetectionMatrix}).
#'
#' @param object A collection of signatures inheriting from "\code{\link{GeneSetCollection}}" or "\code{\link{tbl_geneset}}".
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#' @param method Prediction method. See section "Prediction methods".
#' @param ... Additional arguments affecting the predictions produced.
#'
#' @rdname INTERNAL_predictAnyGeneSetClass
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

    # Update the Hancock metadata.
    # Add global new (general) metadata before existing (method-specific) metadata
    existingHancockMetadata <- metadata(se)[[getPackageName()]]
    newHancockMetadata <- list(
        GeneSets=object,
        method=method,
        packageVersion=packageVersion(getPackageName())
    )
    metadata(se)[[getPackageName()]] <- append(newHancockMetadata, existingHancockMetadata)

    se
}

# predictByProportionPositive ----

#' Identify the Dominant Signatures in Clusters of Samples
#'
#' The \code{predictByProportionPositive} function computes the proportion of samples positive for each signature in each (predefined) cluster
#' and identifies the predominant signature in each cluster.
#' The function stores information tracing the prediction process in the \code{metadata} slot. See Details.
#'
#' @details
#' The function populates the \code{"Hancock"} element of the \code{metadata} slot with the following values:
#' \describe{
#' \item{\code{GeneSets}}{Signatures used to make the predictions}
#' \item{\code{method}}{Name of the method used to make the predictions}
#' \item{\code{packageVersion}}{\code{Hancock} version used to make the predictions}
#' \item{\code{ProportionPositiveByCluster}}{Matrix indicating the proportion of samples in each cluster that are positive for each signature.}
#' \item{\code{TopSignatureByCluster}}{Named vector indicating the predominant signature for each cluster.}
#' }
#'
#' @param object A collection of signatures inheriting from "\code{\link{GeneSetCollection}}" or "\code{\link{tbl_geneset}}".
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param cluster.col Name of a column in \code{colData(se)} that contains
#' a factor indicating cluster membership for each column (i.e. sample) in \code{se}.
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#' @param threshold Value \emph{above which} the marker is considered detected.
#'
#' @return The object \code{se}, updated as follows:
#' \itemize{
#' \item in the \code{metadata} slot, a \code{"Hancock"} item is added (or updated) with information tracing the prediction process. See Details.
#' \item in the \code{colData} slot, a \code{DataFrame} is nested in a new (or updated) \code{"Hancock"} column.
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
#' @seealso \code{\link{predict.GeneSetCollection}}, \code{\link{predict.tbl_geneset}}
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
#' gsc <- GeneSetCollection(list(
#'     GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
#'     GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
#' ))
#' colData(se)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se), replace=TRUE))
#'
#' # Example usage ----
#' library(circlize)
#' # Identify the dominant signature in each cluster
#' se <- predictByProportionPositive(gsc, se, cluster.col="cluster")
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
    uniqueMarkersIds <- uniqueMarkers(object)
    stopifnot(all(uniqueMarkersIds %in% rownames(se)))
    markerDetectionMatrix <- makeMarkerDetectionMatrix(se, uniqueMarkersIds, threshold, assay.type)
    signatureMatrix <- makeSignatureDetectionMatrix(markerDetectionMatrix, object)

    clusterNames <- levels(clusterData)
    numberCellsInCluster <- table(clusterData)

    proportionPositiveByCluster <- matrix(
        data=NA_real_,
        nrow=length(clusterNames),
        ncol=ncol(signatureMatrix),
        dimnames=list(cluster=clusterNames, signature=colnames(signatureMatrix)))
    for (signatureName in colnames(signatureMatrix)) {
        countSignatureInCluster <- table(clusterData, signatureMatrix[, signatureName])[, "TRUE"]
        proportionPositiveByCluster[, signatureName] <- countSignatureInCluster / numberCellsInCluster
    }

    # For each cluster, identify the most frequent signature
    # TODO: warning if ties
    maxColIdx <- max.col(proportionPositiveByCluster, ties.method="first")
    maxSignatureName <- factor(colnames(proportionPositiveByCluster)[maxColIdx])
    names(maxSignatureName) <- rownames(proportionPositiveByCluster)

    # Assign most frequent signature to every cell in each cluster
    newHancockColData <- DataFrame(
        prediction=maxSignatureName[colData(se)[, cluster.col, drop=TRUE]]
    )
    colData(se)[[getPackageName()]] <- newHancockColData

    # Store the proportion of cluster positive for each signature in metadata for later plotting
    newHancockMetadata <- list(
        ProportionPositiveByCluster=proportionPositiveByCluster,
        TopSignatureByCluster=maxSignatureName
    )
    metadata(se)[[getPackageName()]] <- newHancockMetadata

    se
}
