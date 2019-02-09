# Contributing

## Considerate contributions are very welcome!

Here are a few guidelines to help develop and maintain a consistent coding style.

## Table of contents

* [Primary objectives](#primary-objectives)
* [Proof of concepts](#proof-of-concepts)
* [Coding style](#coding-style)
* [Unit tests and code coverage](#unit-tests-and-code-coverage)
* [Internal functions](#internal-functions)
* [New prediction methods](#new-prediction-methods)
* [New plotting functions](#new-plotting-functions)
* [New learning methods](#new-learning-methods)
* [Terminology](#terminology)

## Primary objectives

Most importantly, contributions should directly contribute to the primary objectives of the package, namely:

1. Apply signatures to assign cell identities to new data sets
2. Learn new signatures from data sets, in a format compatible with (1.)

Existing methods implemented in other _R_ packages are welcome.
Those should be handled as described in the [New prediction methods](#new-prediction-methods) section and unit tested as described in the [Unit tests and code coverage](#unit-tests-and-code-coverage) section.

## Proof of concepts

Ideally, a proof-of-concept Rmarkdown notebook should demonstrate the method _before_ adding any new function in `hancock` (i.e., explicitly declaring any new function in the notebook itself).

This can save significant time through community feedback and suggestions from both expert developers and prospective users on the implementation and expected usage of a new functionality _before_ investing significant time and effort into packaging and documenting functions.

For an example, please refer to the proof-of-concept of the function `predictByProportionPositive` available [here](https://github.com/kevinrue/hancock2018/blob/3e065aca67071338b1fcf496790da239ace5425c/1-proportion_signature.Rmd).

Proof-of-concept vignettes may be subsequently updated to demonstrate the same use case, but calling functions implemented in the package. Refer to the demonstration of the function `learnMarkersByPositiveProportionDifference` available [here](https://github.com/kevinrue/hancock2018/blob/f08ee1d34c6bea722757870a339ad3940a48040c/2-learn-signatures.Rmd).


## Coding style

This package adheres to the _Bioconductor_ coding style (https://bioconductor.org/developers/how-to/coding-style/).

Please use common _Bioconductor_ methods and classes, in particular `SummarizedExperiment` and `GeneSetCollection`.
Note that new data structures for gene sets and signatures are under active development.
Those include:
- `BaseSets`  ([_unisets_](https://github.com/kevinrue/unisets) package)
- `tbl_geneset` ([_GeneSet_](https://github.com/Kayla-Morrell/GeneSet) package)

More details are available at https://bioconductor.org/developers/how-to/commonMethodsAndClasses/.

## Unit tests and code coverage

Code coverage should remain at 100%.
Every function, both internal and exported, should be accompanied with its own unit test(s) as part of the _same_ pull request.

A single unit test may include multiple `expect_*` assertions. Use as many `expect_*` as appropriate.

Note that large functions that include several `if` statements and require multiple unit tests to cover every scenario can generally be refactored in multiple smaller functions easier to unit test individually.

The following code is useful to track down lines that are not covered by unit tests:

```
library(covr)
pc <- package_coverage()
report(pc)
```

## NEWS file

Until the package is made available through the [_Bioconductor_](https://bioconductor.org) project, all new features should be described under the "hancock 0.99.0" section, as part of the _same_ pull request.
This will produce a manifest of functionality to accompany the initial submission to the _Bioconductor_ project.

## Internal functions

Internal functions should also be documented using roxygen comments (http://r-pkgs.had.co.nz/man.html) using Markdown formatting (https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html).
However, those do not have to be as comprehensive as exported functions.
Required sections are:

- `@title`
- `@description`
- `@param`
- `@return`
- `@rdname INTERNAL_<...>` with `<...>` being the name of the function (without any trailing ".").
    Make sure that your `.gitignore` contains the entry `INTERNAL_*`. Do _not_ push INTERNAL documentation online.
- `@author`

## New prediction methods

New prediction methods should be first implemented as a separate functions, individually exported in the `NAMESPACE` file.
All prediction methods must accept `object` and `se` as their first two arguments, respectively:

1. the `GeneSetCollection`, `BaseSets`, or `tbl_geneset`
2. the `SummarizedExperiment`

Additional, method-specific parameters may be accepted from the third argument onward.

Once implemented as its own function, a new method should be made available through the [`.predictAnyGeneSetClass`](https://github.com/kevinrue/hancock/blob/c70f7a86026e2ec1d6d54b1d6d732d2a07320c89/R/predict-methods.R#L114) function using a unique `method` identifier.
Make sure the new identifier and method are documented in the `?predictSignatures` man page.

Prediction methods should return the input `SummarizedExperiment` object updated as follows:

- In the `colData` slot, a `DataFrame` nested in a new (or updated) `"hancock"` column should contain at least a first column called `prediction`.
    This column should be populated with the highest-confidence prediction for each sample.
    Additional, method-specific columns may be present from the second column onward.
- In the `metadata` slot, a `list` in a new (or updated) `"hancock"` element, should contain at least the following elements:
    - `"GeneSets"`: the object of class `GeneSetCollection`, `BaseSets`, or `tbl_geneset` containing the signatures used to make the predictions.
    - `"method"`: Identifier of the method used to make the predictions
    - `"packageVersion"`: Version of the `hancock` package used to make the predictions
    - Additional, method-specific elements may appear _after_ the above general metadata

For an example template, please refer to the prediction method [`predictByProportionPositive`](https://github.com/kevinrue/hancock/blob/e7e7f18fb82f59240078533de0c42545485acf9b/R/predict-methods.R#L176), made available using the [`"ProportionPositive"` or `"PP"`](https://github.com/kevinrue/hancock/blob/e7e7f18fb82f59240078533de0c42545485acf9b/R/predict-methods.R#L95) identifiers.

## New plotting functions

New plotting functions should accept `se` as their first argument, namely a `SummarizedExperiment` returned by any prediction method (see above).

Most importantly, plotting function should first check that the input `se` object contains the results of the associated prediction method(s).
In the future, this requirement may be deprecated by the definition of `SummarizedExperiment` subclasses, "flagging" the presence of specific prediction results.

Plotting functions should return a minimal `ggplot2::ggplot` or `ComplexHeatmap::Heatmap` object, giving users maximal freedom to customize the plot.

For an example, please refer to [`plotProportionPositive`](https://github.com/kevinrue/hancock/blob/e7e7f18fb82f59240078533de0c42545485acf9b/R/plot-methods.R#L11), using the result of the [`"ProportionPositive"` or `"PP"`](https://github.com/kevinrue/hancock/blob/e7e7f18fb82f59240078533de0c42545485acf9b/R/predict-methods.R#L176) method.

## New learning methods

Similarly to [new prediction methods](#new-prediction-methods), new learning methods should be first implemented as a separate functions, individually exported in the `NAMESPACE` file.
All prediction methods must accept `se` as their first argument, namely the `SummarizedExperiment` from which to learn signatures.
Additional method-specific parameters may be accepted from the second argument onward.

Once implemented as its own function, a new method should be made available through the [`learnSignatures`](https://github.com/kevinrue/hancock/blob/c70f7a86026e2ec1d6d54b1d6d732d2a07320c89/R/learn-methods.R#L49) method using a unique `method` identifier.
Make sure the new identifier and method are documented in the `?learnSignatures` man page.

Learning methods should return an object inheriting from the `BaseSets` class, defined in the [_unisets_](https://github.com/kevinrue/unisets) package.

For an example template, please refer to the prediction method [`learnMarkersByPositiveProportionDifference`](https://github.com/kevinrue/hancock/blob/master/R/learn-methods.R#L114), made available using the [`"PositiveProportionDifference"` or `"PPD"`](https://github.com/kevinrue/hancock/blob/master/R/learn-methods.R#L54) identifiers.

Metadata produced by learning methods may be stored as [additional columns](https://github.com/kevinrue/hancock/blob/master/R/learn-methods.R#L188) of the `BaseSets` returned.

## Terminology

Until the [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl) or [Human Cell Atlas](https://www.humancellatlas.org) come up with some reference terminology, avoid the use of "cell type" and "(sub-)<sub>n</sub>types" in the code and accompanying documentation.
Those terms are increasingly confusing and open for interpretation as single-cell technologies continuously advance our understanding of _cell differentiation_ into functionally distinct _cell populations_ or _compartments_ currently discriminated by their respective canonical set of _cell surface proteins_ and _transcriptional profiles_ (a few examples of terms that more specifically address individual aspects of the definition of "cell types").
