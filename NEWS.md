# hancock 0.99.0

* Added a `NEWS.md` file to track changes to the package.
* Added Travis continuous integration.
* Added `testthat` unit tests.
* Added `codecov` code coverage.
* Added vignette discussing concepts related to the package functionality.
* Added methods for the detection of markers and signatures:
    `positiveForMarker` (generic), `makeMarkerDetectionMatrix`,
    `makeSignatureDetectionMatrix`.
* Added support for gene set classes
    `GeneSetCollection`, `tbl_geneset`, and `BaseSets`
    defined in packages `GSEABase`, `GeneSet`, and `unisets`.
* Added prediction method `ProportionPositive`.
* Added learning method `PositiveProportionDifference`.
* Added plotting function `plotProportionPositive`.
* Added helper functions: `makeFilterExpression`, `uniqueMarkers`,
    `uniqueSetNames`.
* Added _Shiny_ function: `shinyLabels`.
