
# Example data ----

ncells <- 100
u <- matrix(rpois(20000, 2), ncol=ncells)
rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))
se <- SummarizedExperiment(assays=list(counts=u))

# learnSignatures ----

test_that("learnSignatures works for method ProportionPositive", {
    dummyCluster <- factor(sample(head(LETTERS, 3), ncol(se), replace=TRUE))
    colData(se)[, "cluster"] <- dummyCluster
    nMarkersPerCluster <- 2L

    out <- learnSignatures(se, method="PositiveProportionDifference", cluster.col="cluster", n=nMarkersPerCluster)

    expect_s3_class(out, "tbl_geneset")
    expect_identical(nrow(out), nMarkersPerCluster*nlevels(dummyCluster))
})

test_that("learnPositiveMarkersByProportionDifference requires argument col.cluster", {
    expect_error(
        learnPositiveMarkersByProportionDifference(se),
        "cluster.col is required for method 'PositiveProportionDifference'"
    )
})

test_that("learnPositiveMarkersByProportionDifference requires detection rate in the range 0-1", {
    expect_error(
        learnPositiveMarkersByProportionDifference(se, cluster.col="cluster", min.diff=10),
        "Detection rates are computed in the range 0-1."
    )
})
