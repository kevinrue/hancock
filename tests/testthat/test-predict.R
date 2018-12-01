
# Example data ----

ncells <- 100
u <- matrix(rpois(20000, 2), ncol=ncells)
rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))
se <- SummarizedExperiment(assays=list(counts=u))

gsc <- GeneSetCollection(list(
    GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
    GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
))

# predict.GeneSetCollection ----

test_that("predict.GeneSetCollection works for method ProportionPositive", {
    dummyCluster <- factor(sample(head(LETTERS, 3), ncol(se), replace = TRUE))
    colData(se)[, "cluster"] <- dummyCluster

    out <- predict(gsc, se, method="ProportionPositive", cluster.col="cluster")

    expect_s4_class(out$Hancock, "DataFrame")
    expect_named(out$Hancock, c("prediction"))
    expect_s3_class(out$Hancock$prediction, "factor")
    expect_true(all(out$Hancock$prediction %in% names(gsc)))

    expect_identical(
        names(metadata(out)[["Hancock"]]),
        c("GeneSetCollection", "method", "packageVersion", "ProportionPositiveByCluster",  "TopSignatureByCluster")
    )

    expect_s4_class(metadata(out)[["Hancock"]][["GeneSetCollection"]], "GeneSetCollection")

    expect_identical(metadata(out)[["Hancock"]][["method"]], "ProportionPositive")

    expect_s3_class(metadata(out)[["Hancock"]][["packageVersion"]], "package_version")

    ProportionPositiveByCluster <- metadata(out)[["Hancock"]][["ProportionPositiveByCluster"]]
    expect_is(ProportionPositiveByCluster, "matrix")
    expect_identical(nrow(ProportionPositiveByCluster), nlevels(se$cluster))
    expect_identical(ncol(ProportionPositiveByCluster), length(gsc))

    TopSignatureByCluster <- metadata(out)[["Hancock"]][["TopSignatureByCluster"]]
    expect_s3_class(TopSignatureByCluster, "factor")
    expect_length(TopSignatureByCluster, nlevels(se$cluster))
})

test_that("predictProportionSignatureByCluster requires argument col.cluster", {
    expect_error(
        predictProportionSignatureByCluster(gsc, se),
        "cluster.col is required for method 'ProportionPositive'"
    )
})
