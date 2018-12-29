
# Example data ----

tgs <- tbl_geneset(
    "Cell type 1"=c("Gene001", "Gene002"),
    "Cell type 2"=c("Gene002", "Gene003", "Gene004")
)

# .makeFilterExpressionFromGeneSetCollection ----

test_that(".makeFilterExpressionFromGeneSetCollection works", {
    out <- makeFilterExpression(tgs)

    expect_named(out, c("Cell type 1", "Cell type 2"))
})

# uniqueMarkers ----

test_that("uniqueMarkers works", {
    out <- uniqueMarkers(tgs)
    expect_identical(out, c("Gene001", "Gene002", "Gene003", "Gene004"))
})

# uniqueSetNames ----

test_that("uniqueSetNames works", {
    out <- uniqueSetNames(tgs)
    expect_identical(out, c("Cell type 1", "Cell type 2"))
})
