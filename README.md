
# AIPW

<!-- badges: start -->
<!-- badges: end -->

Augmented inverse probability weighting (AIPW) for binary exposure and outcome

## Installation

``` r
devtools::install_github("yqzhong7/AIPW")
```

## Example

The current version of AIPW was built on the influence function from TMLE,

``` r
library(AIPW)
library(SuperLearner)
library(TMLE)

a <- read.csv("./data.csv",stringsAsFactors = F)
glmnet_learner = create.Learner("SL.glmnet",tune=list(alpha = seq(0,1)))
sl.lib <- c(glmnet_learner$names,"SL.mean","SL.glm")

#tmle
W = subset(a,select=c(high_school,married,employed,white,age,BMI,BPS))
fit <- tmle(Y=as.numeric(a$outcome=="live birth"),
            A=a$treatment,
            W=data.frame(W),
            Q.SL.library=sl.lib,
            g.SL.library=sl.lib,
            family="binomial")

#exposure and outcome must be binary numeric vectors
exposure <- a$treatment
outcome <- as.numeric(a$outcome=="live birth")

aipw(tmle_fit = fit,exposure=exposure,outcome=outcome)
```

