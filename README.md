
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


``` r
library(AIPW)
library(SuperLearner)

#setup data
N <- 200
outcome <- rbinom(N,1,0.3)
exposure <- rbinom(N,1,0.5)
covariates.Q <- matrix(c(rbinom(N,1,0.4),
                                     rnorm(N,mean = 0,sd=1),
                                     rpois(N,lambda = 2)),
                                   ncol=3)

covariates.g <- matrix(c(rbinom(N,1,0.4),
                         rnorm(N,mean = 0,sd=1),
                         rpois(N,lambda = 2)),
                       ncol=3)
                                     
# covariates.g <- c(rbinom(N,1,0.4)) #a vector of a single covariate is also supported
              
#SuperLearner libraries for outcome and exposure models              
sl.Q.lib <-sl.g.lib <- sl.lib <- c("SL.mean","SL.glm")
              
#get individual estimates               
aipw_input_value<- aipw_input(Y=outcome,
                              A=exposure,
                              W.Q=covariates.Q,
                              W.g=covariates.g,
                              Q.SL.library=sl.Q.lib,
                              g.SL.library=sl.g.lib,
                              k_split = 2,
                              verbose = T)
#estimate average treatment effect                                 
aipw(aipw_input = aipw_input_value)                  
```

Use TMLE fitted object as input (sample splitting & separte sets of covariates are not supported using TMLE fitted object):

```R
library(tmle)
tmle object
fit <- tmle(Y=outcome,
            A=exposure,
            W=covariates.Q,
            Q.SL.library=sl.lib,
            g.SL.library=sl.lib,
            family="binomial")

#exposure and outcome must be binary numeric vectors
aipw(tmle_fit = fit, A=exposure, Y=outcome)
```

