STOCfree: prediction of probabilities of freedom from infection from
longitudinal data
================

The aim of the STOCfree package is to predict herd level probabilities
of freedom from infection from longitudinal data collected as part of
surveillance programmes.

Below, the use of the package is presented through the analysis of a
small dataset collected for the surveillance of infection by the BVDV
virus in cattle.

# Package installation and update

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

In order to install (or update) the STOCfree package, run the following
line:

``` r
install_github("AurMad/STOCfree")
```

# Attaching packages

The `STOCfree` package needs to be attached.

``` r
library(STOCfree)
```

The list of available functions and datasets can be accessed by typing

``` r
help(package="STOCfree")
```

We also attach the following packages that will be used later:

``` r
library(ggplot2)
```

# Data

The `STOCfree` package contains a dataset called `herdBTM` which
contains the results of antibody ELISA tests performed on bulk tank
milk. Each row is a testing date in a herd. There are 100 herds with 11
tests for each herd.

There are also data for 2 risk factors of new infection.

  - `ln_nOrig6_12` is the number of source herds for cattle purchased
    between 6 and 12 months before the test
  - `LocalSeroPrev` is the municipality level prevalence of seropositive
    herds on the previous round of testing

<!-- end list -->

``` r
head(herdBTM)
```

    ## # A tibble: 6 x 6
    ##   Farm  DateOfTest   ODR TestResult ln_nOrig6_12 LocalSeroPrev
    ##   <chr> <date>     <dbl>      <int>        <dbl>         <dbl>
    ## 1 FR001 2011-09-20  79.2          1            0          0.2 
    ## 2 FR001 2012-01-12  70.5          1            0          0.1 
    ## 3 FR001 2012-09-25  63.9          1            0          0.2 
    ## 4 FR001 2013-02-05  54.2          0            0          0.18
    ## 5 FR001 2013-08-29  75.8          1            0          0.18
    ## 6 FR001 2014-02-04  67.1          1            0          0.12

# The STOCfree\_data class

For running the model the data needs to be put in a special format: the
STOCfree\_data class. The aim is to help setting the different
parameters required and to help checking the different hypotheses that
are made. The creation of a STOCfree\_data is realised using the
`STOCfree_data()` function.

## Test and risk factor data

In order to demonstrate how to include categorical data in the model, we
create a categorical variable called from `purch_yn` (for purchase
yes/no) from the `ln_nOrig6_12` variable.

``` r
data(herdBTM)
herdBTM$purch_yn <- ifelse(herdBTM$ln_nOrig6_12 > 0, 1, 0)
herdBTM
```

    ## # A tibble: 1,100 x 7
    ##    Farm  DateOfTest   ODR TestResult ln_nOrig6_12 LocalSeroPrev purch_yn
    ##    <chr> <date>     <dbl>      <int>        <dbl>         <dbl>    <dbl>
    ##  1 FR001 2011-09-20  79.2          1        0              0.2         0
    ##  2 FR001 2012-01-12  70.5          1        0              0.1         0
    ##  3 FR001 2012-09-25  63.9          1        0              0.2         0
    ##  4 FR001 2013-02-05  54.2          0        0              0.18        0
    ##  5 FR001 2013-08-29  75.8          1        0              0.18        0
    ##  6 FR001 2014-02-04  67.1          1        0              0.12        0
    ##  7 FR001 2014-09-10  54.4          0        0              0.1         0
    ##  8 FR001 2015-02-01  47.1          0        0              0.08        0
    ##  9 FR001 2015-09-03  39.8          0        0              0.1         0
    ## 10 FR001 2016-02-18  40.5          0        0.693          0.22        1
    ## # … with 1,090 more rows

The `herdBTM` dataset is put in the `STOCfree_data` using the
`STOCfree_data()`. The function takes the following arguments:

  - `test_data`: a data.frame with test results
  - `test_herd_col`: name of the column containing the herd identifier
  - `test_date_col`: name of the column containing the date of test.
    Should be a data formatted as YYYY-mm-dd.
  - `test_res_col`: name of the column containing test results
  - `test_level`: level of testing. Can be ‘herd’ or ‘animal’ (will be
    implemented soon)
  - `risk_factor_data`: a data.frame containing the risk factor for each
    test. Can be the same dataset as above.
  - `risk_herd_col`: name of the column containing the herd identifier
  - `risk_date_col`: name of the date containing the dates on which risk
    factors act on the probability of new infection
  - `risk_factor_col`: name of the column(s) containing the risk factors
  - `risk_factor_type`: can be ‘categorical’ or ‘continuous’. Should be
    a vector of the same length as ‘risk\_factor\_col’

Below, a dataset with the first 3 herds of the herdBTM dataset is put in
the STOCfree data format.

``` r
sfd <- STOCfree_data(test_data = herdBTM[herdBTM$Farm %in% c("FR001", "FR002", "FR003"),],
                          test_herd_col = "Farm",
                          test_date_col = "DateOfTest",
                          test_res_col = "TestResult",
                          test_level = "herd",
                          risk_factor_data = herdBTM[herdBTM$Farm %in% c("FR001", "FR002", "FR003"),],
                          risk_herd_col = "Farm",
                          risk_date_col = "DateOfTest",
                          risk_factor_col = c("purch_yn", "LocalSeroPrev"),
                          risk_factor_type = c("categorical", "continuous"))
```

## Priors for test characteristics

The model accounts for imperfect tests i.e. with test sensitivity and
specificity that are below one. Hypotheses about the test
characteristics is incorporated in the form of priors with Beta
distributions. The current value for the parameters of these Beta
distributions can be checked with the `show_tests()` function.

``` r
show_tests(sfd)
```

    ##         test Se_a Se_b Sp_a Sp_b
    ## 1 TestResult   NA   NA   NA   NA

The values are currently NA and need to be set to reasonable values.

``` r
sfd <- set_priors_tests(sfd,
                 Se_a = 12,
                 Se_b = 2,
                 Sp_a = 200,
                 Sp_b = 4)
```

We can check that the STOCfree\_data object has been updated.

``` r
show_tests(sfd)
```

    ##         test Se_a Se_b Sp_a Sp_b
    ## 1 TestResult   12    2  200    4

The `plot_priors_tests()` function is used to plot the priors.

``` r
plot_priors_tests(sfd)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Priors for infection dynamics

Priors on 2 parameters related to infection dynamics need to be set in
order for the model to work. These are the probability of infection on
the first test and the probability of eliminating the infection for
herds that were infected on the previous month. Both priors are
represented with Beta distributions.

Again, there is a `show_inf_dyn()` that allows displaying the current
values for the prior distributions.

``` r
show_inf_dyn(sfd)
```

    ##  pi1_a  pi1_b tau2_a tau2_b 
    ##     NA     NA     NA     NA

The priors are set:

``` r
sfd <- set_priors_inf_dyn(sfd,
                   pi1_a = 1,
                   pi1_b = 2,
                   tau2_a = 30, 
                   tau2_b = 2)
```

The prior distributions are plotted.

``` r
plot_priors_inf_dyn(sfd)
```

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Priors for risk factors of new infection

From the information provided when constructin the STOCfree\_data
object, the following risk factors have been listed:

``` r
show_rf(sfd)
```

    ##      risk_factor        type modality ref n_obs mean_prior sd_prior
    ## 1      Intercept   intercept       NA   0    33         NA       NA
    ## 3  LocalSeroPrev  continuous       NA   0    33         NA       NA
    ## 11      purch_yn categorical        0   1    29         NA       NA
    ## 2       purch_yn categorical        1   0     4         NA       NA

The association between these risk factors and the probability of new
infection between consecutive month is represented by normal
distributions on a logit scale. The mean and standard deviation of these
distributions need to be set using the `set_priors_rf()` function. The
function needs to be used for each parameter in turn.

``` r
sfd <- set_priors_rf(sfd,
                   risk_factor = "Intercept",
                   mean = -2, sd = 1)

sfd <- set_priors_rf(sfd,
                   risk_factor = "LocalSeroPrev",
                   mean = 0, sd = 2)

sfd <- set_priors_rf(sfd,
                   risk_factor = "purch_yn",
                   modality = 1,
                   mean = 0, sd = 2)

show_rf(sfd)
```

    ##      risk_factor        type modality ref n_obs mean_prior sd_prior
    ## 1      Intercept   intercept       NA   0    33         -2        1
    ## 3  LocalSeroPrev  continuous       NA   0    33          0        2
    ## 11      purch_yn categorical        0   1    29         NA       NA
    ## 2       purch_yn categorical        1   0     4          0        2

The distributions for these priors can then be plotted using the
`plot_priors_rf()` function.

``` r
plot_priors_rf(sfd)
```

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

The STOCfree\_data object is now ready to be modelled.

# STOC free model

## Compiling

The model is run in JAGS. Before extracting samples from the parameter
posterior distributions, it is first compiled. The function takes a
STOCfree\_data object as first argument as well as the number of chains
wanted.

``` r
compiled_model <- compile_JAGS(sfd, n_chains = 4)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 32
    ##    Unobserved stochastic nodes: 196
    ##    Total graph size: 1914
    ## 
    ## Initializing model

## Samples from the parameter posterior distributions

Samples from the posterior distributions of model parameters, predicted
probabilities of infection and predicted statuses are drawn using the
`sample_model` function.

``` r
samples <- sample_model(compiled_model, 
                        n_burnin = 100, 
                        n_iter = 100, 
                        n_thin = 5)
```

# Results

The model retruns a list with 2 components:

  - samples from the model parameters posterior distributions
  - samples from the predicted probability of infection posterior
    distributions

## Model parameters

The samples from the model parameter posterior distributions are in the
`parameters` component of the variable in which we have stored the model
results. We put them in a new variable called `param` to explore these
distributions.

``` r
param <- samples$parameters

param
```

    ## # A tibble: 80 x 9
    ##    .chain .iteration .draw    Se    Sp  tau2 theta.1 theta.2 theta.3
    ##     <int>      <int> <int> <dbl> <dbl> <dbl>   <dbl>   <dbl>   <dbl>
    ##  1      1          1     1 0.854 0.984 0.861   -4.39   3.86   -2.82 
    ##  2      1          2     2 0.854 0.988 0.908   -3.59   0.534   0.876
    ##  3      1          3     3 0.805 0.982 0.934   -4.06  -0.320  -3.62 
    ##  4      1          4     4 0.940 0.978 0.968   -4.50  -3.96   -0.239
    ##  5      1          5     5 0.837 0.981 0.961   -4.61  -2.49   -2.12 
    ##  6      1          6     6 0.882 0.976 0.932   -3.31   1.32   -1.17 
    ##  7      1          7     7 0.849 0.988 0.949   -4.53   3.63   -2.11 
    ##  8      1          8     8 0.780 0.957 0.942   -3.32  -1.13   -0.233
    ##  9      1          9     9 0.956 0.972 0.963   -4.46  -1.62   -1.33 
    ## 10      1         10    10 0.764 0.967 0.949   -4.57   2.18    0.274
    ## # … with 70 more rows

The columns of this dataset are:

  - `.chain`: chain number. Between 1 and 4 in the example because we
    ran 4 MCMC chains in JAGS. This number is the value we selected for
    `n_chains` in the call to `compile_JAGS`.
  - `.iteration`: iteration number. Defined by the n\_iter value in the
    call to sample\_model
  - `.draw`: number of iterations from the first iteration in chain 1,
    pooling all the chains
  - `Se`: sample for sensitivity for this draw
  - `Sp`: sample for pecificity for this draw
  - `tau2`: sample for tau2 for this draw
  - `theta.1`: sample for the logistic model intercept for this draw
  - `theta.2`: sample for the coefficient associated with the first risk
    factor for this draw
  - `theta.3`: sample for the coefficient associated with the second
    risk factor for this draw

![](README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## Probability of infection

The samples from posterior distributions of the predicted probabilities
of infection are in the `proba_inf` component of the variable in which
we have stored the model results. We put them in a new variable called
`proba_inf` to explore these distributions.

``` r
proba_inf <- samples$proba_inf

proba_inf
```

    ## # A tibble: 240 x 6
    ## # Groups:   herd_id [3]
    ##    .chain .iteration .draw herd_id predicted_proba predicted_status
    ##     <int>      <int> <int>   <int>           <dbl>            <dbl>
    ##  1      1          1     1       1         0.0122                 0
    ##  2      1          1     1       2         0.00313                0
    ##  3      1          1     1       3         0.0122                 0
    ##  4      1          2     2       1         0.0268                 0
    ##  5      1          2     2       2         0.00437                0
    ##  6      1          2     2       3         0.0268                 0
    ##  7      1          3     3       1         0.0169                 0
    ##  8      1          3     3       2         0.00325                0
    ##  9      1          3     3       3         0.0169                 0
    ## 10      1          4     4       1         0.0110                 0
    ## # … with 230 more rows

The first columns of this dataset are the same as above. The following
columns differ:

  - `herd_id`: herd identifier
  - `predicted_proba`: predicted probability of infection on this draw
  - `predicted_status`: predicted status regarding infection on this
    draw. This status was simulated by using `predicted_proba` for this
    draw as the proportion parameter in a Bernoulli distribution.

Below we represent the densities for the predicted probabilities of
infection for both the draws predicted as infected and the draws
predicted as uninfected. From this, a rule to categorise herds as free
from infection or not needs to be created.

``` r
cols <- hcl(h = c(120, 0), c = 100, l = 85)

ggplot(proba_inf, aes(x = predicted_proba, colour = factor(predicted_status))) +
  geom_density() +
  xlim(0, 1) +
  xlab("Predicted probability of infection") +
  scale_color_manual(labels = c("Uninfected", "Infected"), values = cols) +
  guides(colour = guide_legend(title="Predicted status"))
```

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Below are the probability densities for each of the herds modelled. Each
color represents a herd.

``` r
ggplot(proba_inf, aes(x = predicted_proba, colour = factor(herd_id))) +
  geom_density() +
  xlim(0, 1) +
  xlab("Predicted probability of infection") +
  theme(legend.position="none")
```

![](README_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

The steps envisonned to categorise herds as free from infection are:

  - define a summary for each density, e.g. 5th percentile of the
    predicted probability of infection
  - choose a cut-off to define freedom from infection, e.g. categorise
    herds with the 5th percentile below 0.1 as free from infection
