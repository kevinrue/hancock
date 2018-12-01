
gsc <- GeneSetCollection(list(
    GeneSet(setName="Cell type 1", c("Gene001", "Gene002")),
    GeneSet(setName="Cell type 2", c("Gene003", "Gene004"))
))

test_that(".makeFilterExpressionFromGeneSetCollection works", {

    out <- .makeFilterExpressionFromGeneSetCollection(gsc)

    expect_named(out, c("Cell type 1", "Cell type 2"))

})
