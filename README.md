AIPW: Augmented Inverse Probability Weighting
================

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/yqzhong7/AIPW/branch/master/graph/badge.svg)](https://codecov.io/gh/yqzhong7/AIPW?branch=master)
[![Travis build
status](https://travis-ci.com/yqzhong7/AIPW.svg?branch=master)](https://travis-ci.com/yqzhong7/AIPW)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#WIP)
[![R build
status](https://github.com/yqzhong7/AIPW/workflows/R-CMD-check/badge.svg)](https://github.com/yqzhong7/AIPW/actions)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

**Authors:** [Yongqi Zhong](https://github.com/yqzhong7), [Ashley
Naimi](https://github.com/ainaimi), [Gabriel
Conzuelo](https://github.com/gconzuelo)

## Installation

``` r
install.packages("remotes")
remotes::install_github("yqzhong7/AIPW")
```

-----

## Example

### Setup example data

``` r
set.seed(888)
N <- 200
outcome <- rbinom(N,1,0.3)
exposure <- rbinom(N,1,0.2)
#covaraites for outcome model (Q)
covariates.Q <- matrix(c(rbinom(N,1,0.4),
                                rnorm(N,mean = 0,sd=1),
                                rpois(N,lambda = 2)),
                              ncol=3)
#covariates for exposure model (g)
covariates.g <- matrix(c(rbinom(N,1,0.4),
                                rnorm(N,mean = 0,sd=1),
                                rpois(N,lambda = 2)),
                              ncol=3)

# covariates.g <- c(rbinom(N,1,0.4)) #a vector of a single covariate is also supported
```

### Create an AIPW object

#### Use [SuperLearner](https://cran.r-project.org/web/packages/SuperLearner/index.html) libraries (reference: [Guide to SuperLearner](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html))

``` r
library(AIPW)
library(SuperLearner)

#SuperLearner libraries for outcome (Q) and exposure models (g)
sl.lib <- c("SL.mean","SL.glm")

#construct an aipw object for later estimations 
AIPW_SL <- AIPW$new(Y= outcome,
                    A= exposure,
                    W.Q=covariates.Q, 
                    W.g=covariates.g,
                    Q.SL.library = sl.lib,
                    g.SL.library = sl.lib,
                    k_split = 3,
                    verbose=FALSE)
```

#### Use [sl3](https://tlverse.org/sl3/index.html) libraries (reference: [Intro to sl3](https://tlverse.org/sl3/articles/intro_sl3.html))

``` r
library(AIPW)
library(sl3)

##construct sl3 learners for outcome (Q) and exposure models (g)
lrnr_glm <- Lrnr_glm$new()
lrnr_mean <- Lrnr_mean$new()
#stacking two learner (this will yield estimates for each learner)
stacklearner <- Stack$new(lrnr_glm, lrnr_mean) 
#metalearner is required to combine the estimates from stacklearner
metalearner <- Lrnr_nnls$new()
sl3.lib <- Lrnr_sl$new(learners = stacklearner,
                      metalearner = metalearner)

#construct an aipw object for later estimations 
AIPW_sl3 <- AIPW$new(Y= outcome,
                    A= exposure,
                    W.Q=covariates.Q, 
                    W.g=covariates.g,
                    Q.SL.library = sl3.lib,
                    g.SL.library = sl3.lib,
                    k_split = 3,
                    verbose=FALSE)
```

### Fit the data to the AIPW object

``` r
#fit the AIPW_SL object
AIPW_SL$fit()
# #fit the AIPW_sl3 object
# AIPW_sl3$fit()
```

### Calculate average treatment effects

#### Check the balance of propensity scores by exposure status

``` r
library(ggplot2)
AIPW_SL$plot.p_score()
#> ATE has not been calculated.
```

![](man/figures/ps_raw-1.png)<!-- -->

#### Estimate the ATE with propensity scores truncation

``` r
#estimate the average causal effects from the fitted AIPW_SL object 
AIPW_SL$calculate_result(g.bound = 0.25) #propensity score truncation 
#>                 Estimate    SE 95% LCL  95% UCL   N
#> Risk Difference   -0.267 0.132 -0.5250 -0.00817 200
#> Risk Ratio         0.393 0.413  0.1748  0.88221 200
#> Odds Ratio         0.266 0.616  0.0797  0.88957 200

#estimate the average causal effects from the fitted AIPW_sl3 object 
# AIPW_sl3$calculate_result(g.bound = 0.25) #propensity score truncation 
```

#### Check the balance of propensity scores by exposure status after truncation

``` r
AIPW_SL$plot.p_score()
```

![](man/figures/ps_trunc-1.png)<!-- -->
