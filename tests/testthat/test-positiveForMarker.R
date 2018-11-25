
# Setup data for this set of tests ----

ncells <- 100
u <- matrix(rpois(20000, 2), ncol=ncells)
rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))

# positiveForMarker ----

test_that("positiveForMarker works for matrix", {
    
    out <- positiveForMarker(u, "Gene001", 0)
    
    expected <- (u["Gene001", ] > 0)
    expect_identical(out, expected)
})

test_that("positiveForMarker works for matrix", {
    
    sce <- SummarizedExperiment(assays=list(counts=u))
    
    out <- positiveForMarker(sce, "Gene001", 0, assay="counts")
    
    expected <- (assay(sce, "counts")["Gene001", ] > 0)
    expect_identical(out, expected)
    
    # Catch invalid assay names
    expect_error(
        positiveForMarker(sce, "Gene001", 0, assay="test")
    )
    
})
