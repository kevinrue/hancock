# Hancock 0.99.0

* Added a `NEWS.md` file to track changes to the package.
* Added Travis CI.
* Added `testthat` unit tests.
* Added `codecov` code coverage.
* Added vignette discussing concepts.
* Added methods for the detection of markers and signatures:
    `positiveForMarker` (generic), `makeMarkerDetectionMatrix`,
    `makeSignatureDetectionMatrix`.
* Added methods for the application of signatures:
    `predict.GeneSetCollection`, `predictProportionSignatureByCluster`.
* Added predicting method `ProportionPositive`.
* Added plotting method `plotProportionPositive`.
* Added various internal helper functions.
