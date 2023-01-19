
# snSMART

<!-- badges: start -->

[![License: GPL
v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![R-CMD-check](https://github.com/sidiwang/osbma/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sidiwang/osbma/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The aim of the **osbma** R package is to predict overall survival in
solid tumor oncology studies using joint models combined through
Bayesian model averaging.

In solid tumor oncology studies, overall survival (OS) is the definitive
endpoint for accessing patient benefit. Reliable OS predictions are very
valuable for both drug development and patients’ end-of-life medical
care. Existing methods often ignore the underlying process of tumor
burden measurements and incur “information loss”. Motivated by an
advanced renal cell carcinoma (RCC) clinical trial, in this paper, we
propose a multivariate joint modeling approach to assess the underlying
dynamics of the progression-free survival (PFS) components to forecast
the death times of trial participants. Through Bayesian model averaging,
our proposed method improves the accuracy of the OS forecast by
combining joint models developed based on each granular component of
PFS. A case study of the RCC trial is conducted, and our method provides
the most accurate predictions across all tested scenarios. The
reliability of our proposed method is also verified through extensive
simulation studies, which include the scenario where OS is completely
independent of PFS. Overall, the proposed methodology provides a
promising candidate for reliable OS prediction in solid tumor oncology
studies.

## Installation

Install the development version from GitHub:

``` r
# Install devtools first if you haven't done so
library(devtools)
# install osbma
devtools::install_github("sidiwang/osbma")
library(osbma)
```

## **osbma** functions covered in this package

- The main function that generates posterior samples based on user input
  trial dataset - `osbma`
- Overall survival prediction based on posterior samples generated
  through MCMC - `predict.osbma`
- Plots generated based on predicted overall survival -
  `plot.predict.data`

## Example

`osbma`: We call the `osbma` function using data from an oncology trial
with 400 total individuals.

``` r
data <- osbma::data
result <- osbma::osbma(data = data, covariate = "trt", method = "rjags")
```

    ## Loading required namespace: rjags

    ## Compiling rjags model...
    ## Calling the simulation using the rjags method...
    ## Adapting the model for 8000 iterations...
    ## Burning in the model for 1000 iterations...
    ## Running the model for 10000 iterations...
    ## Simulation complete
    ## Note: Summary statistics were not produced as there are >50 monitored
    ## variables
    ## [To override this behaviour see ?add.summary and ?runjags.options]
    ## FALSEFinished running the simulation
    ## Calculating summary statistics...
    ## Calculating the Gelman-Rubin statistic for 1629 variables....
    ## Compiling rjags model...
    ## Calling the simulation using the rjags method...
    ## Adapting the model for 8000 iterations...
    ## Burning in the model for 1000 iterations...
    ## Running the model for 10000 iterations...
    ## Simulation complete
    ## Note: Summary statistics were not produced as there are >50 monitored
    ## variables
    ## [To override this behaviour see ?add.summary and ?runjags.options]
    ## FALSEFinished running the simulation

``` r
options(max.print = 40)
summary(result)
```

    ## 
    ## Iterations = 9001:18991
    ## Thinning interval = 10 
    ## Number of chains = 2 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##                         Mean        SD  Naive SE Time-series SE
    ## beta_0[1]          4.425e+01 9.995e-01 2.235e-02      1.090e-01
    ## beta_0[2]         -3.238e-01 8.461e-02 1.892e-03      7.784e-03
    ## beta_1[1]          4.658e+00 1.242e-01 2.778e-03      1.351e-02
    ## beta_1[2]         -1.171e+00 7.172e-02 1.604e-03      7.647e-03
    ## lambda             5.669e-03 6.753e-04 1.510e-05      3.477e-05
    ## alpha_os.1         2.690e+00 2.525e-01 5.646e-03      1.460e-02
    ## eta_2              2.113e-01 1.151e-04 2.573e-06      0.000e+00
    ## alpha_2            1.279e-01 2.242e-02 5.013e-04      0.000e+00
    ## beta_2[1]         -2.806e-03 2.596e-03 5.805e-05      0.000e+00
    ## beta_2[2]         -1.320e-02 2.609e-03 5.834e-05      0.000e+00
    ##  [ reached getOption("max.print") -- omitted 2023 rows ]
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##                         2.5%         25%         50%         75%       97.5%
    ## beta_0[1]          4.236e+01   4.353e+01   4.423e+01   4.495e+01   4.621e+01
    ## beta_0[2]         -4.914e-01  -3.796e-01  -3.242e-01  -2.680e-01  -1.496e-01
    ## beta_1[1]          4.386e+00   4.576e+00   4.661e+00   4.747e+00   4.881e+00
    ## beta_1[2]         -1.303e+00  -1.221e+00  -1.173e+00  -1.126e+00  -1.022e+00
    ## lambda             4.401e-03   5.198e-03   5.652e-03   6.135e-03   6.971e-03
    ## alpha_os.1         2.274e+00   2.513e+00   2.668e+00   2.846e+00   3.213e+00
    ## eta_2              2.112e-01   2.112e-01   2.113e-01   2.114e-01   2.114e-01
    ## alpha_2            1.055e-01   1.055e-01   1.279e-01   1.504e-01   1.504e-01
    ##  [ reached getOption("max.print") -- omitted 2025 rows ]

`predict.osbma`: Here, we call the prediction function to predict the
overall survival of all patients:

``` r
prediction <- predict(result)
```

`plot.predict.osbma`: The plot function provides 3 types of plots: - the
Kaplan-Meier curve of overall survival

``` r
plot(prediction, trt.col.name = "trt", type = "KM")
```

    ## Warning: `select_()` was deprecated in dplyr 0.7.0.
    ## ℹ Please use `select()` instead.
    ## ℹ The deprecated feature was likely used in the dplyr package.
    ##   Please report the issue at <]8;;https://github.com/tidyverse/dplyr/issueshttps://github.com/tidyverse/dplyr/issues]8;;>.

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> - the
patients’ survival timeline of the trial

``` r
# plot(prediction, trt.col.name = "trt", type = "date")
```

- the estimated time of the nth death in the trial (here we use the
  120th death as an example)

``` r
# plot(prediction, trt.col.name = "trt", type = 120)
```
