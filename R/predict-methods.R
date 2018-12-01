
#' Predict Method for GeneSetCollection signatures
#'
#' @param object A set of signatures of class inheriting from "\code{\link{GeneSetCollection}}".
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#' @param method Prediction method. See section "Prediction methods".
#' @param ... Additional arguments affecting the predictions produced.
#'
#' @section Prediction methods:
#' Cluster-level predictions:
#' \describe{
#' \item{ProportionPositive}{
#' Computes the proportion cells positive for each signature in each cluster.
#' Assigns to each cluster the signature detected in the highest proportion of cells.}
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
#'
#' @importFrom S4Vectors metadata
#'
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' ncells <- 100
#' u <- matrix(rpois(20000, 2), ncol=ncells)
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
#' colData(se1)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se1), replace = TRUE))
#' predict(gsc, se1, method="ProportionPositive", cluster.col="cluster")
predict.GeneSetCollection <- function(
    object, se, assay.type="counts", method=c("ProportionPositive"), ...
) {
    method <- match.arg(method)

    if (identical(method, "ProportionPositive")) {
        se <- predictProportionPositiveForGeneSetCollection(object, se, assay.type, ...)
    }

    # Update the Hancock metadata.
    # Add global new (general) metadata before existing (method-specific) metadata
    existingHancockMetadata <- metadata(se)[[getPackageName()]]
    newHancockMetadata <- list(
        GeneSetCollection=object,
        method=method,
        packageVersion=packageVersion(getPackageName())
    )
    metadata(se)[[getPackageName()]] <- append(newHancockMetadata, existingHancockMetadata)

    se
}

#' @param cluster.col Name of column in \code{colData(se)} that contains a factor indicating cluster membership for each column (i.e. sample) in \code{se}.
#' @param threshold Value \emph{above which} the marker is considered detected.
#'
#' @rdname predict.GeneSetCollection
#'
#' @importFrom SummarizedExperiment colData colData<- metadata<-
#' @importFrom S4Vectors DataFrame
#'
#' @author Kevin Rue-Albrecht
predictProportionPositiveForGeneSetCollection <- function(
    object, se, assay.type, cluster.col, threshold=0
) {
    # Sanity checks
    stopifnot(!missing(cluster.col))
    clusterData <- colData(se)[, cluster.col, drop=TRUE]
    stopifnot(is.factor(clusterData))

    # Compute the proportion of each cluster positive for each signature
    uniqueMarkers <- .uniqueMarkers(object)
    stopifnot(all(uniqueMarkers %in% rownames(se)))
    markerDetectionMatrix <- makeMarkerDetectionMatrix(se, uniqueMarkers, threshold, assay.type)
    signatureMatrix <- makeSignatureDetectionMatrix(markerDetectionMatrix, object)

    clusterNames <- levels(clusterData)
    numberCellsInCluster <- table(clusterData)

    proportionPositiveByCluster <- matrix(
        data=FALSE,
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

#' Extract Unique Markers From a GeneSetCollection
#'
#' @rdname INTERNAL_uniqueMarkers
#'
#' @param object A GeneSetCollection object.
#' @author Kevin Rue-Albrecht
.uniqueMarkers <- function(
    object
) {
    # NOTE: later, we may trim signatures to features present in `se`
    # NOTE: in which case, signatures trimmed to length 0 would have to be dropped (!)
    uniqueMarkers <- unique(unlist(geneIds(object)))
    uniqueMarkers
}


#' Identify Markers and Signatures Present in Individual Samples
#'
#' @rdname makeDetectionMatrices
#'
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param markers A character vector, subset of
#' @param threshold Value \emph{above which} the marker is considered detected.
#' @param assay.type A string specifying which assay values to use, e.g., "\code{counts}" or "\code{logcounts}".
#'
#' @return
#' \describe{
#' \item{\code{makeMarkerDetectionMatrix}}{A logical matrix indicating the presence of each marker in each sample.}
#' \item{\code{makeSignatureDetectionMatrix}}{A logical matrix indicating the presence of each signature in each sample.}
#' }
#'
#' @export
#'
#' @author Kevin Rue-Albrecht
#'
#' @examples
#' # Example data ----
#' library(SummarizedExperiment)
#' ncells <- 100
#' u <- matrix(rpois(20000, 2), ncol=ncells)
#' rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
#' se <- SummarizedExperiment(assays=list(counts=u))
#'
#' gsc <- GeneSetCollection(list(
#'   GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
#'   GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
#' ))
#'
#' # Example usage ----
#' markerMatrix <- makeMarkerDetectionMatrix(se, unique(unlist(geneIds(gsc))))
#' signatureMatrix <- makeSignatureDetectionMatrix(markerMatrix, gsc)
makeMarkerDetectionMatrix <- function(
    se, markers, threshold=0, assay.type="counts"
) {
    if (any(duplicated(markers))) {
        warning("Dropping duplicated markers values")
        markers <- unique(markers)
    }
    markerDetectionMatrix <- matrix(
        data=FALSE,
        nrow=ncol(se), ncol=length(markers),
        dimnames=list(colnames(se), markers))
    # TODO: lapply, BiocParallel, ...
    for (gene in markers) {
        markerDetectionMatrix[, gene] <- positiveForMarker(se, gene, threshold, assay.type)
    }
    markerDetectionMatrix
}

#' @rdname makeDetectionMatrices
#'
#' @param matrix A logical matrix indicating the presence of each marker in each sample.
#' See \code{\link{makeMarkerDetectionMatrix}}
#' @param object A set of signatures of class inheriting from "\code{\link{GeneSetCollection}}".
#'
#' @export
#'
#' @importFrom S4Vectors FilterRules evalSeparately
makeSignatureDetectionMatrix <- function(matrix, object) {
    filterExpressions <- .makeFilterExpressionFromGeneSetCollection(object)
    fr <- FilterRules(filterExpressions)
    es <- evalSeparately(fr, as.data.frame(matrix))
    es
}
