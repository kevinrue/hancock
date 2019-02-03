
# Example data ----

geneLists <- list(
    "Cell type 1"=c("Gene001", "Gene002"),
    "Cell type 2"=c("Gene002", "Gene003", "Gene004")
)

bs <- as(geneLists, "BaseSets")

# .makeFilterExpressionFromGeneSetCollection ----

test_that(".makeFilterExpressionFromGeneSetCollection works", {
    out <- makeFilterExpression(bs)

    expect_named(out, c("Cell type 1", "Cell type 2"))
})

# uniqueMarkers ----

test_that("uniqueMarkers works", {
    out <- uniqueMarkerNames(bs)
    expect_identical(out, c("Gene001", "Gene002", "Gene003", "Gene004"))
})

# uniqueSetNames ----

test_that("uniqueSetNames works", {
    out <- uniqueSetNames(bs)
    expect_identical(out, c("Cell type 1", "Cell type 2"))
})
