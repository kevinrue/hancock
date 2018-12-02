
## Considerate contributions are more than welcome!

Here are a few guidelines to help develop and maintain a consistent coding style.

## Primary objectives

Most importantly, contributions should directly contribute to the primary objectives of the package, namely:

1. Apply signatures to assign cell identities to new data sets
2. Learn new signatures from data sets, in a format compatible with (1.)

## Proof of concepts

Ideally, a proof-of-concept Rmarkdown notebook should demonstrate the method _before_ adding any new function in `Hancock` (i.e., explicitly declaring any new function in the notebook itself).

This can save significant time through community feedback and suggestions from both expert developers and prospective users on the implementation and expected usage of a new functionality _before_ investing significant time and effort into packaging and documenting functions.

For an example, please refer to the proof-of-concept of the `predictProportionSignatureByCluster` function available [here](https://github.com/kevinrue/Hancock2018/blob/master/1-proportion_signature.Rmd).


## Coding style

This package follows the _Bioconductor_ coding style (https://bioconductor.org/developers/how-to/coding-style/).

Use common _Bioconductor_ methods and classes, in particular `SummarizedExperiment` and `GeneSetCollection`.
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
All prediction methods must accept `object` and `se` as their first two arguments, respectively the `GeneSetCollection` and `SummarizedExperiment` used to make predictions.
Additional method-specific parameters may be accepted from the third argument onward.

Once implemented as its own function, a new method should be made available through the `predict.GeneSetCollection` function using a unique `method` identifier.

Prediction methods should return the input `SummarizedExperiment` object updated as follows:

- In the `colData` slot, a `DataFrame` nested in a new (or updated) `"Hancock"` column should contain at least a first column called `prediction`. Additional, method-specific columns may be present from the second column onward.
- In the `metadata` slot, a `list` in a new (or updated) `"Hancock"` element, should contain at least the following elements:
    - `"GeneSetCollection"`: the `GeneSetCollection` object used to make the predictions
    - `"method"`: Identifier of the method used to make the predictions
    - `"packageVersion"`: Version of the `Hancock` package used to make the predictions
    - Additional, method-specific elements may appear _after_ the above general metadata

For an example template, please refer to the prediction method `predictProportionSignatureByCluster`, made available using the `"ProportionPositive"` identifier.

## New plotting functions

New plotting functions should accept `se` as their first argument, namely a `SummarizedExperiment` returned by any prediction method (see above).

Most importantly, plotting function should first check that the input `se` object contains the results of the associated prediction method(s).

Plotting functions should return a minimal `ggplot2::ggplot` or `ComplexHeatmap::Heatmap` object, giving users maximal freedom to customize the plot.

For an example, please refer to `plotProportionPositive`, using the result of the `predictProportionSignatureByCluster` method.
