<img src="inst/www/hancock_hexsticker.png" align="right" alt="" width="120" />

[![Travis build status](https://travis-ci.org/kevinrue/hancock.svg?branch=master)](https://travis-ci.org/kevinrue/hancock)
[![Coverage status](https://codecov.io/gh/kevinrue/hancock/branch/master/graph/badge.svg)](https://codecov.io/github/kevinrue/hancock?branch=master)

# hancock

The goal of the [_hancock_](https://github.com/kevinrue/hancock) package is to provide a collection of methods for learning and applying gene signatures associated with cellular phenotypes and identities.
Particular focus is given to single-cell data stored in objects derived from the [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) class.

# Prerequisites

The [_hancock_](https://github.com/kevinrue/hancock) package supports classes of gene sets defined in multiple packages.
However, it uses classes defined in the `unisets` to return newly learned signatures with accompanying metadata.
This dependency may be installed as follows:

```
install.packages("devtools")
devtools::install_github("kevinrue/unisets", build_opts = c("--no-resave-data", "--no-manual"))
```

Several functions support the `tbl_geneset` class defined in the `GeneSet` package.
This package is currently only available from [GitHub](https://github.com/Kayla-Morrell/GeneSet).
It may be installed as follows:

```
devtools::install_github("Kayla-Morrell/GeneSet")
```

# Installation

The [_hancock_](https://github.com/kevinrue/hancock) package may be installed as follows:

```
install.packages("devtools")
devtools::install_github("kevinrue/hancock")
```

To install the vignette as well (building it requires an additional minute or so), please use the following code:

```
devtools::install_github("kevinrue/hancock", build_opts = c("--no-resave-data", "--no-manual"))
```

# Usage

Demonstration notebooks are available on the companion repository: https://github.com/kevinrue/hancock2018

# Contributing

Considerate contributions are very welcome!
Please refer to the [contributing guidelines](https://github.com/kevinrue/hancock/blob/master/CONTRIBUTING.md) for more details.
