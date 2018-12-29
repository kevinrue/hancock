
tgs <- tbl_geneset(
    "Cell type 1"=c("Gene001", "Gene002"),
    "Cell type 2"=c("Gene002", "Gene003", "Gene004")
)

# .panelGeneration ----

test_that(".panelGeneration requires predictByProportionPositive results", {
    out <- .panelGeneration(tgs)
    expect_s3_class(out, "shiny.tag")
})
