AIPW: Augmented Inverse Probability Weighting
================

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/yqzhong7/AIPW/branch/master/graph/badge.svg)](https://codecov.io/gh/yqzhong7/AIPW?branch=master)
[![Travis build
status](https://travis-ci.com/yqzhong7/AIPW.svg?branch=master)](https://travis-ci.com/yqzhong7/AIPW)
[![R build
status](https://github.com/yqzhong7/AIPW/workflows/R-CMD-check/badge.svg)](https://github.com/yqzhong7/AIPW/actions)
[![](https://www.r-pkg.org/badges/version/AIPW?color=blue)](https://cran.r-project.org/package=AIPW)

<!-- badges: end -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

**Contributors:** [Yongqi Zhong](https://github.com/yqzhong7), [Ashley
Naimi](https://github.com/ainaimi), [Gabriel
Conzuelo](https://github.com/gconzuelo), [Edward
Kennedy](https://github.com/ehkennedy)

------------------------------------------------------------------------

Augmented inverse probability weighting (AIPW) is a doubly robust
estimator for causal inference. The `AIPW` package is designed for
estimating the average treatment effect of a binary exposure on risk
difference (RD), risk ratio (RR) and odds ratio (OR) scales with
user-defined stacked machine learning algorithms
([SuperLearner](https://CRAN.R-project.org/package=SuperLearner) or
[sl3](https://tlverse.org/sl3/index.html)). Users need to examine causal
assumptions (e.g., consistency) before using this package.

If you find this package is helpful, please consider to cite:

    @article{zhong_aipw_2021,
        author = {Zhong, Yongqi and Kennedy, Edward H and Bodnar, Lisa M and Naimi, Ashley I},
        title = {AIPW: An R Package for Augmented Inverse Probability Weighted Estimation of Average Causal Effects},
        journal = {American Journal of Epidemiology},
        year = {2021},
        month = {07},
        issn = {0002-9262},
        doi = {10.1093/aje/kwab207},
        url = {https://doi.org/10.1093/aje/kwab207},
    }

------------------------------------------------------------------------

## Contents:

-   ##### [Installation](#Installation)

-   ##### [Example](#Example)

    -   ###### [Setup example data](#data)

    -   ###### [One line version](#one_line)

-   ##### [Parallelization and progress bar](#par)

-   ##### [Use tmle/tmle3 as input](#tmle)

-   ##### [References](#ref)

------------------------------------------------------------------------

## <a id="Installation"></a>Installation

### CRAN version

``` r
install.packages("AIPW")
```

### Github version

``` r
install.packages("remotes")
remotes::install_github("yqzhong7/AIPW")
```

**\* CRAN version only supports SuperLearner and tmle. Please install
the Github version (master branch) if you choose to use sl3 and tmle3.**

## <a id="Example"></a>Example

### <a id="data"></a>Setup example data

``` r
set.seed(888)
data("eager_sim_obs")
outcome <- eager_sim_obs$sim_Y
exposure <- eager_sim_obs$sim_A
#covariates for both outcome model (Q) and exposure model (g)
covariates <- as.matrix(eager_sim_obs[-1:-2])

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
AIPW_SL <- AIPW$new(Y = outcome,
                    A = exposure,
                    W = covariates, 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 3,
                    verbose=FALSE)$
  fit()$
  #Default truncation is set to 0.025; using 0.25 here is for illustrative purposes and not recommended
  summary(g.bound = 0.25)$ 
  plot.p_score()$
  plot.ip_weights()
```

![](man/figures/one_line-1.png)<!-- -->![](man/figures/one_line-2.png)<!-- -->

To see the results, set `verbose = TRUE`(default) or:

``` r
print(AIPW_SL$result, digits = 2)
#>                  Estimate    SE 95% LCL 95% UCL   N
#> Risk of exposure     0.44 0.046  0.3528    0.53 118
#> Risk of control      0.31 0.051  0.2061    0.41  82
#> Risk Difference      0.14 0.068  0.0048    0.27 200
#> Risk Ratio           1.45 0.191  0.9974    2.11 200
#> Odds Ratio           1.81 0.295  1.0144    3.22 200
```

To obtain average treatment effect among the treated/controls (ATT/ATC),
`statified_fit()` must be used:

``` r
AIPW_SL_att <- AIPW$new(Y = outcome,
                    A = exposure,
                    W = covariates, 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 3,
                    verbose=T)
suppressWarnings({
  AIPW_SL_att$stratified_fit()$summary()
})
#> Done!
#>                     Estimate     SE  95% LCL 95% UCL   N
#> Risk of exposure      0.4352 0.0467  0.34362   0.527 118
#> Risk of control       0.3244 0.0513  0.22385   0.425  82
#> Risk Difference       0.1108 0.0684 -0.02320   0.245 200
#> Risk Ratio            1.3416 0.1858  0.93210   1.931 200
#> Odds Ratio            1.6048 0.2927  0.90429   2.848 200
#> ATT Risk Difference   0.0991 0.0880 -0.07339   0.272 200
#> ATC Risk Difference   0.1148 0.0634 -0.00946   0.239 200
```

You can also use the `aipw_wrapper()` to wrap `new()`, `fit()` and
`summary()` together (also support method chaining):

``` r
AIPW_SL <- aipw_wrapper(Y = outcome,
                        A = exposure,
                        W = covariates, 
                        Q.SL.library = c("SL.mean","SL.glm"),
                        g.SL.library = c("SL.mean","SL.glm"),
                        k_split = 3,
                        verbose=TRUE,
                        stratified_fit=F)$plot.p_score()$plot.ip_weights()
```

## <a id="par"></a>Parallelization with `future.apply` and progress bar with `progressr`

In default setting, the `AIPW$fit()` method will be run sequentially.
The current version of AIPW package supports parallel processing
implemented by
[future.apply](https://github.com/HenrikBengtsson/future.apply) package
under the [future](https://github.com/HenrikBengtsson/future) framework.
Simply use `future::plan()` to enable parallelization and `set.seed()`
to take care of the random number generation (RNG) problem:

``` r
###Additional steps for parallel processing###
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

###Same procedure for AIPW as described above###
AIPW_SL <- AIPW$new(Y = outcome,
                    A = exposure,
                    W = covariates, 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 3,
                    verbose=TRUE)$fit()$summary()
#> Done!
#>                  Estimate     SE 95% LCL 95% UCL   N
#> Risk of exposure    0.443 0.0462 0.35284   0.534 118
#> Risk of control     0.306 0.0510 0.20607   0.406  82
#> Risk Difference     0.137 0.0677 0.00482   0.270 200
#> Risk Ratio          1.449 0.1906 0.99741   2.106 200
#> Odds Ratio          1.807 0.2946 1.01442   3.219 200
```

Progress bar that supports parallel processing is available in the
`AIPW$fit()` method through the API from
[progressr](https://github.com/HenrikBengtsson/progressr) package:

``` r
library(progressr)
#define the type of progress bar
handlers("progress")
#reporting through progressr::with_progress() which is embedded in the AIPW$fit() method
with_progress({
  AIPW_SL <- AIPW$new(Y = outcome,
                    A = exposure,
                    W = covariates, 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 3,
                    verbose=FALSE)$fit()$summary()
})
#also available for the wrapper
with_progress({
  AIPW_SL <- aipw_wrapper(Y = outcome,
                        A = exposure,
                        W = covariates, 
                        Q.SL.library = c("SL.mean","SL.glm"),
                        g.SL.library = c("SL.mean","SL.glm"),
                        k_split = 3,
                        verbose=FALSE)
})
```

## <a id="tmle"></a>Use `tmle`/`tmle3` fitted object as input (`AIPW_tmle` class)

`AIPW_tmle` class is designed for using `tmle`/`tmle3` fitted object as
input

#### 1. `tmle`

``` r
require(tmle)
#> Loading required package: tmle
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loaded glmnet 4.0
#> Welcome to the tmle package, version 1.4.0.1
#> 
#> Use tmleNews() to see details on changes and bug fixes
require(SuperLearner)
tmle_fit <- tmle(Y = outcome, A = exposure,W = covariates,
                 Q.SL.library=c("SL.mean","SL.glm"),
                 g.SL.library=c("SL.mean","SL.glm"),
                 family="binomial")
tmle_fit
#>  Additive Effect
#>    Parameter Estimate:  0.13484
#>    Estimated Variance:  0.004308
#>               p-value:  0.039939
#>     95% Conf Interval: (0.0061944, 0.26348) 
#> 
#>  Additive Effect among the Treated
#>    Parameter Estimate:  0.13525
#>    Estimated Variance:  0.0043779
#>               p-value:  0.040941
#>     95% Conf Interval: (0.0055666, 0.26494) 
#> 
#>  Additive Effect among the Controls
#>    Parameter Estimate:  0.1332
#>    Estimated Variance:  0.0042363
#>               p-value:  0.040711
#>     95% Conf Interval: (0.0056273, 0.26077) 
#> 
#>  Relative Risk
#>    Parameter Estimate:  1.4352
#>               p-value:  0.051898
#>     95% Conf Interval: (0.99703, 2.0658) 
#> 
#>               log(RR):  0.36128
#>     variance(log(RR)):  0.034539 
#> 
#>  Odds Ratio
#>    Parameter Estimate:  1.7837
#>               p-value:  0.045351
#>     95% Conf Interval: (1.012, 3.1436) 
#> 
#>               log(OR):  0.57866
#>     variance(log(OR)):  0.083597
#extract fitted tmle object to AIPW
AIPW_tmle$
  new(A=exposure,Y=outcome,tmle_fit = tmle_fit,verbose = TRUE)$
  summary(g.bound=0.025)
#> Cross-fitting is supported only within the outcome model from a fitted tmle object (with cvQinit = TRUE)
#>                  Estimate     SE 95% LCL 95% UCL   N
#> Risk of exposure    0.445 0.0454 0.35577   0.534 118
#> Risk of control     0.310 0.0497 0.21243   0.407  82
#> Risk Difference     0.135 0.0656 0.00619   0.263 200
#> Risk Ratio          1.435 0.1815 1.00562   2.048 200
#> Odds Ratio          1.784 0.2818 1.02672   3.099 200
```

#### 2. `tmle3`

``` r
library(sl3)
library(tmle3)
node_list <- list(A = "sim_A",Y = "sim_Y",W = colnames(eager_sim_obs)[-1:-2])
or_spec <- tmle_OR(baseline_level = "0",contrast_level = "1")
tmle_task <- or_spec$make_tmle_task(eager_sim_obs,node_list)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_mean <- make_learner(Lrnr_mean)
sl <- Lrnr_sl$new(learners = list(lrnr_glm,lrnr_mean))
learner_list <- list(A = sl, Y = sl)
tmle3_fit <- tmle3(or_spec, data=eager_sim_obs, node_list, learner_list)

# parse tmle3_fit into AIPW_tmle class
AIPW_tmle$
  new(A=eager_sim_obs$sim_A,Y=eager_sim_obs$sim_Y,tmle_fit = tmle3_fit,verbose = TRUE)$
  summary()
```

------------------------------------------------------------------------

### <a id="ref"></a>References:

Robins JM, Rotnitzky A (1995). Semiparametric efficiency in multivariate
regression models with missing data. Journal of the American Statistical
Association.

Chernozhukov V, Chetverikov V, Demirer M, et al (2018). Double/debiased
machine learning for treatment and structural parameters. The
Econometrics Journal.

Kennedy EH, Sjolander A, Small DS (2015). Semiparametric causal
inference in matched cohort studies. Biometrika.

Pearl, J., 2009. Causality. Cambridge university press.
