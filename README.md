STOCfree: prediction of probabilities of freedom from infection from
longitudinal data
================

The aim of the STOCfree package is to predict herd level probabilities
of freedom from infection from longitudinal data collected as part of
surveillance programmes.

Below, the use of the package is presented through the analysis of a
small dataset collected for the surveillance of infection by the BVDV
virus in cattle.

# Package installation

Before installing the package, you need to make sure that JAGS is
installed. This programme can be installed from the following website:
<https://sourceforge.net/projects/mcmc-jags/files/>.

The easiest way to install the `STOCfree` package is from Github. This
requires installing the `devtool` package first.

``` r
install.packages("devtools")
```

Then load the `devtool` package:

``` r
library(devtools)
```

and install the STOCfree package:

``` r
install_github("AurMad/STOCfree")
```

# Attaching the STOC free package

The `STOCfree` package needs to be attached.

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

herdBTM_month
```

    ## # A tibble: 6,224 x 10
    ##    row_id herd_id date__1    month_id Farm  DateOfTest   ODR TestResult
    ##     <int>   <int> <date>        <int> <chr> <date>     <dbl>      <int>
    ##  1      1       1 2011-09-01        2 FR001 2011-09-20  79.2          1
    ##  2      2       1 2011-10-01        3 <NA>  NA          NA           NA
    ##  3      3       1 2011-11-01        4 <NA>  NA          NA           NA
    ##  4      4       1 2011-12-01        5 <NA>  NA          NA           NA
    ##  5      5       1 2012-01-01        6 FR001 2012-01-12  70.5          1
    ##  6      6       1 2012-02-01        7 <NA>  NA          NA           NA
    ##  7      7       1 2012-03-01        8 <NA>  NA          NA           NA
    ##  8      8       1 2012-04-01        9 <NA>  NA          NA           NA
    ##  9      9       1 2012-05-01       10 <NA>  NA          NA           NA
    ## 10     10       1 2012-06-01       11 <NA>  NA          NA           NA
    ## # … with 6,214 more rows, and 2 more variables: ln_nOrig6_12 <dbl>,
    ## #   LocalSeroPrev <dbl>

# Priors

## Test

Prior for tests are stored in a variable called `test_priors`.

``` r
test_priors <- list(
  Se_beta_a = 12,
  Se_beta_b = 2,
  Sp_beta_a = 200,
  Sp_beta_b = 4
)
```

## Infection dynamics

Probability of being infected on the first testing time for a herd and
probability of not eliminating the infection between 2 consecutive
tests.

``` r
infection_priors <- list(
  pi1_beta_a = 1,
  pi1_beta_b = 2,
  tau2_beta_a = 30,
  tau2_beta_b = 2
 )
```

and priors for risk factors

``` r
risk_factor_priors <- list(
  theta_norm_mean = 0,
  theta_norm_sd = .01
)
```

# STOC free model

## Compiling

The JAGS model is compiled using the `compile_JAGS` function. We just
use 3 herds to test the code.

``` r
test <- expand_month(data = herdBTM[herdBTM$Farm %in% c("FR001", "FR002", "FR003"),],
                     herd_colname = Farm,
                     date_colname = DateOfTest,
                     test_res_colname = TestResult)
```

    ## Joining, by = "date__1"Joining, by = "month_id"Joining, by = c("herd_id",
    ## "date__1")

``` r
compiled_model <- compile_JAGS(test_data = test, 
             herd_id = herd_id, 
             row_id = row_id,
             month = month_id,
             risk_factors = c("ln_nOrig6_12", "LocalSeroPrev"),
             test_result = TestResult,
             test_priors = test_priors, 
             infection_priors = infection_priors, 
             risk_factor_priors = risk_factor_priors,
             n_chains = 2)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 32
    ##    Unobserved stochastic nodes: 196
    ##    Total graph size: 1907
    ## 
    ## Initializing model

## Samples from the parameter posterior distributions

Samples from the posterior distributions of model parameters, predicted
probabilities of infection and predicted statuses are drawn using the
`sample_model` function.

``` r
samples <- sample_model(compiled_model, n_burnin = 100, n_iter = 100, n_thin = 5)
```

# Results

The model retruns a list with 3 components:

  - samples from the model parameters posterior distributions
  - samples from the predicted probability of infection posterior
    distributions
  - predicted model performance

## Model parameters

The samples from the model parameter posterior distributions are stored
in the `parameters` component of the variable in which we have stored
the model results

``` r
param <- samples$parameters
```

## Probability of infection

The samples from the model predicted probability of infection posterior
distributions are stored in the `proba_inf` component of the variable in
which we have stored the model results

``` r
proba_inf <- samples$proba_inf
```

## Model performance

Needs a better description

``` r
cutoff <- samples$choice_cutoff
```
