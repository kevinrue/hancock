# Hancock 0.99.0

* Added a `NEWS.md` file to track changes to the package.
* Added Travis continuous integration.
* Added `testthat` unit tests.
* Added `codecov` code coverage.
* Added vignette discussing concepts related to the package functionality.
* Added methods for the detection of markers and signatures:
    `positiveForMarker` (generic), `makeMarkerDetectionMatrix`,
    `makeSignatureDetectionMatrix`.
* Added support for gene set classes: `GeneSetCollection`, `tbl_geneset`.
* Added predicting method `ProportionPositive`.
* Added plotting method `plotProportionPositive`.
* Added helper functions: `makeFilterExpression`, `uniqueMarkers`.
