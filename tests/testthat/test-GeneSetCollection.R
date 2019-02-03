
# Example data ----

gsc <- GeneSetCollection(list(
    GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
    GeneSet(setName="Cell type 2", c("Gene002", "Gene003", "Gene004"))
))

# .makeFilterExpressionFromGeneSetCollection ----

test_that(".makeFilterExpressionFromGeneSetCollection works", {
    out <- makeFilterExpression(gsc)

    expect_named(out, c("Cell type 1", "Cell type 2"))
})

# uniqueMarkers ----

test_that("uniqueMarkers works", {
    out <- uniqueMarkerNames(gsc)
    expect_identical(out, c("Gene001", "Gene002", "Gene003", "Gene004"))
})

# uniqueSetNames ----

test_that("uniqueSetNames works", {
    out <- uniqueSetNames(gsc)
    expect_identical(out, c("Cell type 1", "Cell type 2"))
})
