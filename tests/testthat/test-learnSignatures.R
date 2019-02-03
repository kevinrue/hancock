
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

    expect_s4_class(out, "BaseSets")
    expect_lte(length(out), nMarkersPerCluster*nlevels(dummyCluster))

    out <- learnSignatures(se, method="PositiveProportionDifference", cluster.col="cluster", n=nMarkersPerCluster, diff.method="min")

    expect_s4_class(out, "BaseSets")
    expect_lte(length(out), nMarkersPerCluster*nlevels(dummyCluster))
})

test_that("learnMarkersByPositiveProportionDifference requires argument col.cluster", {
    expect_error(
        learnMarkersByPositiveProportionDifference(se),
        "cluster.col is required for method 'PositiveProportionDifference'"
    )
})

test_that("learnMarkersByPositiveProportionDifference requires detection rate in the range 0-1", {
    expect_error(
        learnMarkersByPositiveProportionDifference(se, cluster.col="cluster", min.diff=10),
        "min.diff must be a scalar in the range [0,1].",
        fixed=TRUE
    )
})
