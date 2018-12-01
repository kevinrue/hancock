
# Example data ----

gsc <- GeneSetCollection(list(
    GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
    GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
))

# .makeFilterExpressionFromGeneSetCollection ----

test_that(".makeFilterExpressionFromGeneSetCollection works", {
    out <- .makeFilterExpressionFromGeneSetCollection(gsc)

    expect_named(out, c("Cell type 1", "Cell type 2"))
})

# .uniqueMarkers ----

test_that(".uniqueMarkers works", {
    out <- .uniqueMarkers(gsc)
    expect_identical(out, c("Gene001", "Gene002", "Gene003", "Gene004"))
})
