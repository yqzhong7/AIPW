AIPW: Augmented Inverse Probability Weighting
================

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#WIP)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/yqzhong7/AIPW/branch/master/graph/badge.svg)](https://codecov.io/gh/yqzhong7/AIPW?branch=master)
[![Travis build
status](https://travis-ci.com/yqzhong7/AIPW.svg?branch=master)](https://travis-ci.com/yqzhong7/AIPW)
[![R build
status](https://github.com/yqzhong7/AIPW/workflows/R-CMD-check/badge.svg)](https://github.com/yqzhong7/AIPW/actions)

<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

**Authors:** [Yongqi Zhong](https://github.com/yqzhong7), [Ashley
Naimi](https://github.com/ainaimi), [Gabriel
Conzuelo](https://github.com/gconzuelo)

## Contents:

  - ##### [Installation](#Installation)

  - ##### [Example](#Example)
    
      - ###### [Setup example data](#data)
    
      - ###### [One line version](#one_line)

  - ##### [Parallelization](#par)

  - ##### [Use tmle/tmle3 as input](#tmle)

-----

## <a id="Installation"></a>Installation

``` r
install.packages("remotes")
remotes::install_github("yqzhong7/AIPW")
```

## <a id="Example"></a>Example

### <a id="data"></a>Setup example data

``` r
set.seed(888)
N <- 200
outcome <- rbinom(N,1,0.3)
exposure <- rbinom(N,1,0.2)
#covariates for both outcome model (Q) and exposure model (g)
covariates <- matrix(c(rbinom(N,1,0.4),
                                rnorm(N,mean = 0,sd=1),
                                rpois(N,lambda = 2)),
                              ncol=3)

# covariates <- c(rbinom(N,1,0.4)) #a vector of a single covariate is also supported
```

### <a id="one_line"></a>One line version (`AIPW` class: method chaining from R6class)

``` r
library(AIPW)
library(SuperLearner)
#> Loading required package: nnls
#> Super Learner
#> Version: 2.0-26
#> Package created on 2019-10-27
library(ggplot2)
library(progressr) #for the progress bar
AIPW_SL <- AIPW$new(Y= outcome,
                    A= exposure,
                    W=covariates, 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 3,
                    verbose=TRUE)$
  fit()$
  calculate_result(g.bound = 0.25)$
  plot.p_score()
#> Done!
#>                 Estimate    SE 95% LCL 95% UCL   N
#> Risk Difference   -0.253 0.132 -0.5126   0.006 200
#> Risk Ratio         0.417 0.394  0.1927   0.903 200
#> Odds Ratio         0.288 0.600  0.0888   0.934 200
```

![](man/figures/one_line-1.png)<!-- -->

## <a id="par"></a>Parallelization with `future.apply`

In default setting, the `AIPW$fit()` method will be run sequentially.
The current version of AIPW package supports parallel processing
implemented by
[future.apply](https://github.com/HenrikBengtsson/future.apply) package
under the [future](https://github.com/HenrikBengtsson/future) framework.
Simply use `future::plan()` to enable parallelization and `set.seed()`
to take care of the random number generation (RNG) problem:

``` r
# install.packages("future.apply")
library(future.apply)
#> Loading required package: future
future::plan(multiprocess, workers=2, gc=T)
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
set.seed(888)
AIPW_SL <- AIPW$new(Y= outcome,
                    A= exposure,
                    W=covariates, 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 3,
                    verbose=TRUE)$fit()$calculate_result()
#> Done!
#>                 Estimate    SE 95% LCL 95% UCL   N
#> Risk Difference   -0.315 0.167 -0.6425  0.0116 200
#> Risk Ratio         0.357 0.479  0.1395  0.9132 200
#> Odds Ratio         0.220 0.759  0.0498  0.9766 200
```

## <a id="tmle"></a>Use `tmle`/`tmle3` fitted object as input (`AIPW_tmle` class)

`AIPW_tmle` class is designed for using `tmle`/`tmle3` fitted object as
input

#### 1\. `tmle`

As shown in the message,
[tmle](https://cran.r-project.org/web/packages/tmle/index.html) alone
does not support sample splitting. In addition, `tmle` does not support
different covariate sets for the exposure and the outcome models,
respectively.

``` r
require(tmle)
#> Loading required package: tmle
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loaded glmnet 3.0-2
#> Welcome to the tmle package, version 1.4.0.1
#> 
#> Use tmleNews() to see details on changes and bug fixes
require(SuperLearner)
tmle_fit <- tmle(Y=outcome,A=exposure,W=covariates,
                 Q.SL.library=c("SL.mean","SL.glm"),
                 g.SL.library=c("SL.mean","SL.glm"),
                 family="binomial")
AIPW_tmle$
  new(A=exposure,Y=outcome,tmle_fit = tmle_fit,verbose = TRUE)$
  calculate_result(g.bound=0.025)
#> Sample splitting was not supported with a fitted tmle object
#>                 Estimate    SE 95% LCL 95% UCL   N
#> Risk Difference   -0.214 0.158 -0.5235  0.0961 200
#> Risk Ratio         0.438 0.522  0.1575  1.2177 200
#> Odds Ratio         0.326 0.752  0.0745  1.4231 200
```

#### 2\. `tmle3`

Similarly, [tmle3](https://github.com/tlverse/tmle3) may not support
different covariate sets for the exposure and the outcome models,
respectively. However, `tmle3` conducts sample splitting and propensity
truncation (0.025) by default.

``` r
require(tmle3)
#> Loading required package: tmle3
require(sl3)
#> Loading required package: sl3
df <- data.frame(exposure,outcome,covariates)
node_list <- list(A = "exposure",Y = "outcome",W = colnames(df)[-1:-2])
or_spec <- tmle_OR(baseline_level = "0",contrast_level = "1")
tmle_task <- or_spec$make_tmle_task(df,node_list)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_mean <- make_learner(Lrnr_mean)
sl <- Lrnr_sl$new(learners = list(lrnr_glm,lrnr_mean))
learner_list <- list(A = sl, Y = sl)
tmle3_fit <- tmle3(or_spec, data=df, node_list, learner_list)

# parse tmle3_fit into AIPW_tmle class
AIPW_tmle$
  new(A=df$exposure,Y=df$outcome,tmle_fit = tmle3_fit,verbose = TRUE)$
  calculate_result()
#> Propensity scores from fitted tmle3 object are by default truncated (0.025)
#>                 Estimate    SE 95% LCL 95% UCL   N
#> Risk Difference   -0.213 0.158 -0.5236  0.0968 200
#> Risk Ratio         0.438 0.524  0.1570  1.2240 200
#> Odds Ratio         0.326 0.755  0.0743  1.4320 200
```
