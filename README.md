
# AIPW

<!-- badges: start -->
<!-- badges: end -->

Augmented inverse probability weighting (AIPW) for binary exposure and outcome

## Installation

``` r
install.packages("remotes")
remotes::install_github("yqzhong7/AIPW")
```

## Example

The current version of AIPW was built on the influence function from TMLE,

``` r
library(AIPW)
library(SuperLearner)
library(tmle)

#setup
N <- 100
outcome <- rbinom(N,1,0.3)
exposure <- rbinom(N,1,0.5)
covariates <- as.data.frame(matrix(c(rbinom(N,1,0.4),
              rnorm(N,mean = 0,sd=1),
              rpois(N,lambda = 2)
              ),ncol=3))
              
#tmle object
sl.lib <- c("SL.mean","SL.glm")
fit <- tmle(Y=outcome,
            A=exposure,
            W=covariates,
            Q.SL.library=sl.lib,
            g.SL.library=sl.lib,
            family="binomial")

#exposure and outcome must be binary numeric vectors
aipw(tmle_fit = fit, exposure=exposure,outcome=outcome)
```

