
## Considerate contributions are more than welcome!

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

Ideally, a proof-of-concept Rmarkdown notebook should demonstrate the method _before_ adding any new function in `Hancock` (i.e., explicitly declaring any new function in the notebook itself).

This can save significant time through community feedback and suggestions from both expert developers and prospective users on the implementation and expected usage of a new functionality _before_ investing significant time and effort into packaging and documenting functions.

For an example, please refer to the proof-of-concept of the function `predictByProportionPositive` available [here](https://github.com/kevinrue/Hancock2018/blob/3e065aca67071338b1fcf496790da239ace5425c/1-proportion_signature.Rmd).

Proof-of-concept vignettes may be subsequently updated to demonstrate the same use case, but calling functions implemented in the package. Refer to the demonstration of the function `learnMarkersByPositiveProportionDifference` available [here](https://github.com/kevinrue/Hancock2018/blob/f08ee1d34c6bea722757870a339ad3940a48040c/2-learn-signatures.Rmd).


## Coding style

This package adheres to the _Bioconductor_ coding style (https://bioconductor.org/developers/how-to/coding-style/).

Please use common _Bioconductor_ methods and classes, in particular `SummarizedExperiment`, `GeneSetCollection`, and `tbl_geneset` ([GeneSet](https://github.com/Kayla-Morrell/GeneSet) package, in development).
Note that the `GeneSetCollection` class is currently the primary class to store gene sets in the Bioconductor release code; that said, the new `tbl_geneset` class in development is expected to provide a more efficient representation of large data sets and a more familiar `tibble`-like framework.

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

Until the package is made available through the [_Bioconductor_](https://bioconductor.org) project, all new features should be described under the "Hancock 0.99.0" section, as part of the _same_ pull request.
This will produce a manifest of functionality to accompany the initial submission to the _Bioconductor_ project.

## Internal functions

Internal functions should also be documented using roxygen comments (http://r-pkgs.had.co.nz/man.html).
However, those do not have to be as comprehensive as exported functions.
Nevertheless, required sections are:

- A title
- `@rdname INTERNAL_<...>` with `<...>` being the name of the function (without any trailing ".").
    Make sure that your `.gitignore` contains the entry `INTERNAL_*`. Do _not_ push INTERNAL documentation online.
- `@param`
- `@return`
- `@author`

## New prediction methods

New prediction methods should be first implemented as a separate functions, individually exported in the `NAMESPACE` file.
All prediction methods must accept `object` and `se` as their first two arguments, respectively the `GeneSetCollection` or `tbl_geneset`, and `SummarizedExperiment` used to make predictions.
Additional method-specific parameters may be accepted from the third argument onward.

Once implemented as its own function, a new method should be made available through the `predict.GeneSetCollection` function using a unique `method` identifier.
Make sure the new identifier and method are documented in the `?predictHancock` man page.

Prediction methods should return the input `SummarizedExperiment` object updated as follows:

- In the `colData` slot, a `DataFrame` nested in a new (or updated) `"Hancock"` column should contain at least a first column called `prediction`. Additional, method-specific columns may be present from the second column onward.
- In the `metadata` slot, a `list` in a new (or updated) `"Hancock"` element, should contain at least the following elements:
    - `"GeneSets"`: the object of class `GeneSetCollection` or `tbl_geneset` containing the signatures used to make the predictions.
    - `"method"`: Identifier of the method used to make the predictions
    - `"packageVersion"`: Version of the `Hancock` package used to make the predictions
    - Additional, method-specific elements may appear _after_ the above general metadata

For an example template, please refer to the prediction method `predictByProportionPositive`, made available using the `"ProportionPositive"` identifier.

## New plotting functions

New plotting functions should accept `se` as their first argument, namely a `SummarizedExperiment` returned by any prediction method (see above).

Most importantly, plotting function should first check that the input `se` object contains the results of the associated prediction method(s).

Plotting functions should return a minimal `ggplot2::ggplot` or `ComplexHeatmap::Heatmap` object, giving users maximal freedom to customize the plot.

For an example, please refer to `plotProportionPositive`, using the result of the `predictProportionSignatureByCluster` method.

## New learning methods

Similarly to [new prediction methods](#new-prediction-methods), new learning methods should be first implemented as a separate functions, individually exported in the `NAMESPACE` file.
All prediction methods must accept `se` as their first argument, namely the `SummarizedExperiment` from which to learn signatures.
Additional method-specific parameters may be accepted from the second argument onward.

Once implemented as its own function, a new method should be made available through the `learnSignatures` function using a unique `method` identifier.
Make sure the new identifier and method are documented in the `?learnHancock` man page.

Learning methods should return a `tbl_geneset` object, defined in the [GeneSet](https://github.com/Kayla-Morrell/GeneSet) package.

For an example template, please refer to the prediction method `learnMarkersByPositiveProportionDifference`, made available using the `"PositiveProportionDifference"` identifier.

## Terminology

Until the [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl) or [Human Cell Atlas](https://www.humancellatlas.org) come up with some reference terminology, avoid the use of "cell type" and "(sub-)<sub>n</sub>types" in the code and accompanying documentation.
Those terms are increasingly confusing and open for interpretation as single-cell technologies advance our understanding of _cell differentiation_ into functionally distinct _cell populations_ or _compartments_ currently discriminated by their respective canonical set of _cell surface proteins_ and _transcriptional profiles_ (a few examples of terms that more specifically address individual aspects of the definition of "cell types").
