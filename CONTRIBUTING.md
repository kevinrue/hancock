
## Considerate contributions are more than welcome!

Here are a few guidelines to help develop and maintain a consistent coding style.

## Primary objectives

Most importantly, contributions should directly contribute to the primary objectives of the package. namely:

1. Apply signatures to assign cell identities to new data sets
2. Learn new signatures from data sets, in a format compatible with (1.)

## Coding style

This package follows the _Bioconductor_ coding style (https://bioconductor.org/developers/how-to/coding-style/).

Use common _Bioconductor_ methods and classes, in particular `SummarizedExperiment` and `GeneSetCollection`.

## Unit tests and code coverage

Code coverage should remain at 100%.
Every function, both internal and exported, should be accompanied with its own unit test(s).

A single unit test may include multiple `expect_*` assertions. Use as many `expect_*` as appropriate.

Note that large functions that include several `if` statements and require multiple unit tests to cover every scenario can generally be refactored in multiple smaller functions easier to unit test individually.

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
