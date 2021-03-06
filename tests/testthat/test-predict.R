
# Example data ----

ncells <- 100
u <- matrix(rpois(20000, 2), ncol=ncells)
rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))
sce <- SingleCellExperiment(assays=list(counts=u))

reducedDim(sce, "PCA") <- matrix(rnorm(2*ncol(sce)), nrow=ncol(sce), ncol=2)

gsc <- GeneSetCollection(list(
    GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
    GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
))

bs <- as(list(
    "Cell type 1"=c("Gene001", "Gene002"),
    "Cell type 2"=c("Gene003", "Gene004")
), "Sets")

# predict.GeneSetCollection ----

test_that("predict.GeneSetCollection works for method ProportionPositive", {
    dummyCluster <- factor(sample(head(LETTERS, 3), ncol(sce), replace=TRUE))
    colData(sce)[, "cluster"] <- dummyCluster

    out <- predict(gsc, sce, method="ProportionPositive", cluster.col="cluster")

    expect_s4_class(out$hancock, "DataFrame")
    expect_named(out$hancock, c("prediction"))
    expect_s3_class(out$hancock$prediction, "factor")
    expect_true(all(out$hancock$prediction %in% names(gsc)))

    expect_identical(
        names(metadata(out)[["hancock"]]),
        c("GeneSets", "method", "packageVersion", "ProportionPositiveByCluster",  "TopSignatureByCluster")
    )

    expect_s4_class(metadata(out)[["hancock"]][["GeneSets"]], "GeneSetCollection")

    expect_identical(metadata(out)[["hancock"]][["method"]], "ProportionPositive")

    expect_s3_class(metadata(out)[["hancock"]][["packageVersion"]], "package_version")

    ProportionPositiveByCluster <- metadata(out)[["hancock"]][["ProportionPositiveByCluster"]]
    expect_is(ProportionPositiveByCluster, "matrix")
    expect_identical(nrow(ProportionPositiveByCluster), length(gsc))
    expect_identical(ncol(ProportionPositiveByCluster), nlevels(sce$cluster))

    TopSignatureByCluster <- metadata(out)[["hancock"]][["TopSignatureByCluster"]]
    expect_s3_class(TopSignatureByCluster, "factor")
    expect_length(TopSignatureByCluster, nlevels(sce$cluster))

    # Test plotting methods
    plotOut <- plotProportionPositive(out)
    expect_s4_class(plotOut, "Heatmap")

    plotOut <- barplotPredictionCount(out, highlight=c("Cell type 1"), labels=FALSE)
    expect_s3_class(plotOut, "ggplot")

    plotOut <- barplotPredictionProportion(out, highlight=c("Cell type 1"), labels=FALSE)
    expect_s3_class(plotOut, "ggplot")

    plotOut <- reducedDimPrediction(out, highlight=c("Cell type 1"), redDimType="PCA")
    expect_s3_class(plotOut, "ggplot")

    plotOut <- hancock:::.plotWrapper(
        out, highlight=c("Cell type 1"), plotType="reducedDimPrediction", redDimType="PCA")
    expect_s3_class(plotOut, "ggplot")
})

# predict.Sets ----

test_that("predict.Sets works for method ProportionPositive", {
    dummyCluster <- factor(sample(head(LETTERS, 3), ncol(sce), replace=TRUE))
    colData(sce)[, "cluster"] <- dummyCluster

    out <- predict(bs, sce, method="ProportionPositive", cluster.col="cluster")

    expect_s4_class(out$hancock, "DataFrame")
    expect_named(out$hancock, c("prediction"))
    expect_s3_class(out$hancock$prediction, "factor")
    expect_true(all(out$hancock$prediction %in% names(gsc)))

    expect_identical(
        names(metadata(out)[["hancock"]]),
        c("GeneSets", "method", "packageVersion", "ProportionPositiveByCluster",  "TopSignatureByCluster")
    )

    expect_s4_class(metadata(out)[["hancock"]][["GeneSets"]], "Sets")

    expect_identical(metadata(out)[["hancock"]][["method"]], "ProportionPositive")

    expect_s3_class(metadata(out)[["hancock"]][["packageVersion"]], "package_version")

    ProportionPositiveByCluster <- metadata(out)[["hancock"]][["ProportionPositiveByCluster"]]
    expect_is(ProportionPositiveByCluster, "matrix")
    expect_identical(nrow(ProportionPositiveByCluster), length(gsc))
    expect_identical(ncol(ProportionPositiveByCluster), nlevels(sce$cluster))

    TopSignatureByCluster <- metadata(out)[["hancock"]][["TopSignatureByCluster"]]
    expect_s3_class(TopSignatureByCluster, "factor")
    expect_length(TopSignatureByCluster, nlevels(sce$cluster))

    # Test plotting methods
    plotOut <- plotProportionPositive(out)
    expect_s4_class(plotOut, "Heatmap")

    plotOut <- barplotPredictionCount(out, highlight=c("Cell type 1"), labels=FALSE)
    expect_s3_class(plotOut, "ggplot")

    plotOut <- barplotPredictionProportion(out, highlight=c("Cell type 1"), labels=FALSE)
    expect_s3_class(plotOut, "ggplot")

    plotOut <- reducedDimPrediction(out, highlight=c("Cell type 1"), redDimType="PCA")
    expect_s3_class(plotOut, "ggplot")

    plotOut <- hancock:::.plotWrapper(
        out, highlight=c("Cell type 1"), plotType="reducedDimPrediction", redDimType="PCA")
    expect_s3_class(plotOut, "ggplot")
})

# predictByProportionPositive ----

test_that("predictByProportionPositive requires argument col.cluster", {
    expect_error(
        predictByProportionPositive(gsc, sce),
        "cluster.col is required for method 'ProportionPositive'"
    )
})

# plotProportionPositive ----

test_that("plotProportionPositive requires predictByProportionPositive results", {
    expect_error(
        plotProportionPositive(sce),
        "Method 'ProportionPositive' was not run yet.")
})
