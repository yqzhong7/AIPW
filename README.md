
# AIPW

<!-- badges: start -->
<!-- badges: end -->

Augmented inverse probability weighting (AIPW) for binary exposure and outcome

## Installation

``` r
install.packages("remotes")
remotes::install_github("yqzhong7/AIPW",ref="sample-splitting-superlearner")
```

## Example


``` r
library(AIPW)
library(SuperLearner)

#setup data
N <- 100
outcome <- rbinom(N,1,0.3)
exposure <- rbinom(N,1,0.5)
covariates <- as.data.frame(matrix(c(rbinom(N,1,0.4),
              rnorm(N,mean = 0,sd=1),
              rpois(N,lambda = 2)
              ),ncol=3))
              
#SuperLearner libraries for outcome and exposure models              
sl.Q.lib <-sl.g.lib <- sl.lib <- c("SL.mean","SL.glm")
              
#get individual estimates               
aipw_input_value<- aipw_input(Y=outcome,
                                  A=exposure,
                                  W=covariates,
                                  Q.SL.library=sl.Q.lib, #outcome model
                                  g.SL.library=sl.g.lib, #exposure model
                                  k_split = 5)
#estimate average treatment effect                                 
aipw(aipw_input = aipw_input_value)                  
```

Use TMLE fitted object as input (sample splitting is not supported using TMLE fitted object):

```R
library(tmle)
#tmle object
fit <- tmle(Y=outcome,
            A=exposure,
            W=covariates,
            Q.SL.library=sl.lib,
            g.SL.library=sl.lib,
            family="binomial")

#exposure and outcome must be binary numeric vectors
aipw(tmle_fit = fit, A=exposure, Y=outcome)
```

