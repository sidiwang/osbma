
- [osbma](#osbma)
  - [Installation](#installation)
  - [**osbma** functions covered in this
    package](#osbma-functions-covered-in-this-package)
  - [Example](#example)

# osbma

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
    ## Adapting the model for 20000 iterations...
    ## Burning in the model for 10000 iterations...
    ## Running the model for 90000 iterations...
    ## Simulation complete
    ## Calculating summary statistics...
    ## Calculating the Gelman-Rubin statistic for 29 variables....
    ## Finished running the simulation
    ## Compiling rjags model...
    ## Calling the simulation using the rjags method...
    ## Adapting the model for 20000 iterations...
    ## Burning in the model for 10000 iterations...
    ## Running the model for 90000 iterations...
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
    ## Iterations = 30001:119998
    ## Thinning interval = 3 
    ## Number of chains = 2 
    ## Sample size per chain = 30000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##                         Mean        SD  Naive SE Time-series SE
    ## beta_0[1]          2.022e+01 8.211e-01 3.352e-03      2.991e-02
    ## beta_0[2]         -9.837e-02 9.817e-01 4.008e-03      4.109e-03
    ## beta_1[1]         -2.474e+00 2.353e-01 9.605e-04      3.783e-03
    ## beta_1[2]          1.144e+00 1.332e-01 5.440e-04      2.161e-03
    ## lambda            -4.788e-03 1.333e-03 5.442e-06      9.047e-06
    ## alpha_os.1         2.095e+00 1.025e-01 4.183e-04      1.242e-03
    ## eta_2              5.256e-03 5.219e-03 2.131e-05      2.975e-05
    ## alpha_2            5.465e-01 9.007e-02 3.677e-04      4.865e-04
    ## beta_2[1]         -1.520e+00 2.211e-01 9.025e-04      2.775e-03
    ## beta_2[2]          1.382e-02 1.523e-01 6.219e-04      1.913e-03
    ##  [ reached getOption("max.print") -- omitted 3703 rows ]
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##                         2.5%         25%         50%         75%       97.5%
    ## beta_0[1]          1.864e+01   1.966e+01   2.022e+01   2.078e+01   2.182e+01
    ## beta_0[2]         -2.038e+00  -7.579e-01  -9.927e-02   5.667e-01   1.816e+00
    ## beta_1[1]         -2.943e+00  -2.632e+00  -2.471e+00  -2.314e+00  -2.022e+00
    ## beta_1[2]          8.887e-01   1.053e+00   1.142e+00   1.233e+00   1.409e+00
    ## lambda            -7.455e-03  -5.675e-03  -4.773e-03  -3.883e-03  -2.216e-03
    ## alpha_os.1         1.899e+00   2.024e+00   2.094e+00   2.164e+00   2.302e+00
    ## eta_2              1.336e-04   1.518e-03   3.656e-03   7.282e-03   1.936e-02
    ## alpha_2            3.728e-01   4.848e-01   5.450e-01   6.067e-01   7.265e-01
    ##  [ reached getOption("max.print") -- omitted 3705 rows ]

`predict.osbma`: Here, we call the prediction function to predict the
overall survival of all patients:

``` r
prediction <- predict(result)
```

`plot.predict.osbma`: The plot function provides 3 types of plots: 

\* the Kaplan-Meier curve of overall survival

``` r
plot(prediction, trt.col.name = "trt", type = "KM")
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> 

\* the patients’ survival timeline of the trial

``` r
# plot(prediction, trt.col.name = "trt", type = "date")
```

- the estimated time of the nth death in the trial (here we use the
  120th death as an example)

``` r
# plot(prediction, trt.col.name = "trt", type = 120)
```

As `.md` file doesn’t show interactive plots, please check file
`demo.Rmd` for the last two plots. Please download `demo.Rmd` and view
it in browser.
