STOCfree: prediction of probabilities of freedom from infection from
longitudinal data
================

The aim of the STOCfree package is predict herd level probabilities of
freedom from infection from longitudinal data collected as part of
surveillance programmes.

Below, the use of the package is presented through the analysis of a
small dataset collected for the surveillance of infection by the BVDV
virus in cattle.

# Package installation

The easiest way to install the `STOCfree` package is from Github. This
requires installing the `devtool` package first.

``` r
install.packages("devtools")
```

Then load the `devtool` package:

``` r
library(devtools)
```

and install the STOCfree package

``` r
install_github("AurMad/STOCfree")
```

# Loading the STOCfree package

``` r
library(STOCfree)
```

# Data

The `STOCfree` package contains a dataset called `herdBTM` which
contains the results of antibody ELISA tests performed on bulk tank
milk. Each row is a testing date in a herd. There are 100 herds with 11
tests for each herd.

``` r
head(herdBTM)
```

    ##    Farm DateOfTest    ODR TestResult ln_nOrig6_12 LocalSeroPrev
    ## 1 FR001 2011-09-20 79.211          1            0          0.20
    ## 2 FR001 2012-01-12 70.547          1            0          0.10
    ## 3 FR001 2012-09-25 63.907          1            0          0.20
    ## 4 FR001 2013-02-05 54.219          0            0          0.18
    ## 5 FR001 2013-08-29 75.793          1            0          0.18
    ## 6 FR001 2014-02-04 67.068          1            0          0.12

# Formatting the data for analysis

Herds in the `herdBTM` are tested approximately every 6 months. The
STOCfree model models infection with a montlhy time step. The data used
by the model need to have one row per month. The herdBTM data is
expanded to have one row per month with the `expand_month()` function.

``` r
herdBTM_month <- expand_month(data = herdBTM,
                              herd_colname = Farm,
                              date_colname = DateOfTest,
                              test_res_colname = TestResult)
```

# STOC free model

## Compiling

## Burnin

## Sampling

## Results
