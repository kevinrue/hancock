
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

tgs <- tbl_geneset(
    "Cell type 1" = c("Gene001", "Gene002"),
    "Cell type 2" = c("Gene002", "Gene003", "Gene004")
)

# learnSignatures ----

test_that("learnSignatures works for method ProportionPositive", {
    dummyCluster <- factor(sample(head(LETTERS, 3), ncol(se), replace = TRUE))
    colData(se)[, "cluster"] <- dummyCluster
    nMarkersPerCluster <- 2L

    out <- learnSignatures(se, method="ProportionDifference", cluster.col="cluster", n=nMarkersPerCluster)

    expect_s3_class(out, "tbl_geneset")
    expect_identical(nrow(out), nMarkersPerCluster*nlevels(dummyCluster))
})

test_that("learnSignatures requires argument col.cluster", {
    expect_error(
        learnSignaturesByProportionDifference(se),
        "cluster.col is required for method 'ProportionDifference'"
    )
})
