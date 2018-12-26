<img src="inst/www/Hancock_hexsticker.png" align="right" alt="" width="120" />

[![Travis build status](https://travis-ci.org/kevinrue/Hancock.svg?branch=master)](https://travis-ci.org/kevinrue/Hancock)
[![Coverage status](https://codecov.io/gh/kevinrue/Hancock/branch/master/graph/badge.svg)](https://codecov.io/github/kevinrue/Hancock?branch=master)

# Hancock

The goal of Hancock is to provide a collection of methods for learning and applying gene signatures associated with cellular phenotypes and identities.
Particular focus is given to single-cell data stored in objects derived from the [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) class.

# Prerequisites

Several functions depend on the `GeneSet` package, currently only available from [GitHub](https://github.com/Kayla-Morrell/GeneSet).
The dependencies may be installed as follows:

```
install.packages("devtools")
devtools::install_github("Kayla-Morrell/GeneSet", "tibble_implement")
```

# Installation

The `Hancock` package may be installed as follows:

```
install.packages("devtools")
devtools::install_github("kevinrue/Hancock")
```

# Usage

Demonstration notebooks are available on the companion repository: https://github.com/kevinrue/Hancock2018

# Contributing

Considerate contributions are more than welcome.
Please refer to the [contributing guidelines](https://github.com/kevinrue/Hancock/blob/master/CONTRIBUTING.md) for more details.
