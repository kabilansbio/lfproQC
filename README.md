
<!-- README.md is generated from README.Rmd. Please edit that file -->

## About

<!-- badges: start -->
<!-- badges: end -->

The goal of lfproQC R package is to provide an optimal combination of
normalization and imputation methods for the label-free proteomics
expression dataset. Users can also access this R package through the
Shiny application available at <http://omics.icar.gov.in/lfproQC>.

## Installation

You can install the development version of lfproQC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kabilansbio/lfproQC", build_manual = TRUE, build_vignette = TRUE)
```

## Example

This is a basic example for finding the best combinations of
normalization and imputation method for the label-free proteomics
expression dataset:

``` r
library(lfproQC)
## basic example code with the example dataset and data groups
yeast <- best_combination(yeast_data, yeast_groups)
yeast$`Best Combinations`
#> NULL
```

**The overall workflow for using the ‘lfproQC’ package**

<img src="vignettes/images/flow1.jpg" width="800px" />
