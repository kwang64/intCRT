
# BCPP.Package

<!-- badges: start -->
<!-- badges: end -->

## Overview
intCRT R package provides tools for fitting and conducting inference on marginal Cox proportional hazard models for 
cluster-dependent and interval-censored time-to-event data. It is structured according to the method described in the 
Cook et al.(2023), including the Botswana Combination Prevention Project as an applied example.

The package includes two sets of functions:

- functions to estimate the Cox ph model parameters (including both the regression coefficients and the baseline hazard functions) and
estimate their standard error

- functions to generate simulated clustered interval-censored dataset consistent with the structure of
the Botswana Combination Prevention Project

## Installation

You can install the development version of BCPP.Package from [GitHub](https://github.com/) with:

``` r
# Install devtools if necessary
install.packages("devtools")

# Install BCPP from GitHub
devtools::install_github("yourusername/BCPP")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BCPP.Package)
## basic example code
```

