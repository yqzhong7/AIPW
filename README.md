
# AIPW

<!-- badges: start -->
<!-- badges: end -->

Augmented inverse probability weighting (AIPW) for binary exposure and outcome


### Installation

``` r
install.packages("remotes")
remotes::install_github("yqzhong7/AIPW",ref="R6")
```

------

### Example

1. Setup example data

``` r
set.seed(888)
N <- 200
outcome <- rbinom(N,1,0.3)
exposure <- rbinom(N,1,0.5)
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

2. Use [SuperLearner](https://cran.r-project.org/web/packages/SuperLearner/index.html) libraries (reference: [Guide to SuperLearner](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html))

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
                    verbose=T)
#estimate the average causal effects
AIPW_SL$calculate_result()

#example output
#                 Estimate     SE 95% LCL 95% UCL   N
# Risk Difference   -0.150 0.0654 -0.2786 -0.0221 200
# Risk Ratio         0.631 0.6203  0.1872  2.1299 200
# Odds Ratio         0.504 1.3261  0.0374  6.7747 200
```

3. Use [sl3](https://tlverse.org/sl3/index.html) libraries (reference: [Intro to sl3](https://tlverse.org/sl3/articles/intro_sl3.html))


```R
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
                    verbose=T)
#estimate the average causal effects
AIPW_sl3$calculate_result()

#example output
#                 Estimate     SE 95% LCL 95% UCL   N
# Risk Difference   -0.147 0.0659 -0.2763 -0.0179 200
# Risk Ratio         0.638 0.6216  0.1885  2.1559 200
# Odds Ratio         0.511 1.3356  0.0373  7.0037 200
```

