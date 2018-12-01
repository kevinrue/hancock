
# Setup data for this set of tests ----

ncells <- 100
u <- matrix(rpois(20000, 2), ncol=ncells)
rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))

se <- SummarizedExperiment(assays=list(counts=u))

# positiveForMarker ----

test_that("positiveForMarker works for SummarizedExperiment", {
    out <- positiveForMarker(se, "Gene001", 0, assay="counts")

    expected <- (assay(se, "counts")["Gene001", ] > 0)
    expect_identical(out, expected)

    # Catch invalid assay names
    expect_error(
        positiveForMarker(se, "Gene001", 0, assay="test")
    )
})

# makeMarkerDetectionMatrix ----

test_that("makeMarkerDetectionMatrix works", {
    markers <- c("Gene001", "Gene002", "Gene003", "Gene004")
    out <- makeMarkerDetectionMatrix(se, markers, threshold=0, assay.type="counts")

    expect_identical(ncol(out), length(markers))
    expect_identical(nrow(out), ncol(se))
    expect_type(out, "logical")
})

test_that("makeMarkerDetectionMatrix warns about duplicated markers", {

    markers <- c("Gene001", "Gene001")
    expect_warning(
        makeMarkerDetectionMatrix(se, markers, threshold=0, assay.type="counts"),
        "Dropping duplicated markers values"
    )

})

# makeSignatureDetectionMatrix ----

test_that("makeSignatureDetectionMatrix works", {
    nmarkers <- 5
    markerMatrix <- matrix(sample(c(TRUE, FALSE), 100, TRUE), ncol=nmarkers, nrow=ncells)
    colnames(markerMatrix) <- paste0("Marker", sprintf("%02d", seq_len(ncol(markerMatrix))))
    gsc <- GeneSetCollection(
        GeneSet(c("Marker01", "Marker02"), setName="Set1"),
        GeneSet(c("Marker03", "Marker04"), setName="Set2")
    )

    out <- makeSignatureDetectionMatrix(markerMatrix, gsc)
    expect_identical(ncol(out), length(gsc))
    expect_identical(nrow(out), nrow(markerMatrix))
    expect_type(out, "logical")
})
