## INSTALL AND LOAD PACKAGES
packages <- c("data.table","SuperLearner","ggplot2",
              "tidyverse","splines","tmle")
for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos='http://lib.stat.cmu.edu/R/CRAN') 
  }
}

for (package in packages) {
  library(package, character.only=T)
}

## FORMAT FIGURES
thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.title=element_blank(),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

a <- read_csv("./eager_base_imputed.csv")[,-1]
a

xgboost_learner = create.Learner("SL.xgboost",params = list(minobspernode = 20),
                                 tune=list(ntrees = c(500,1000),
                                           max_depth = c(4,5),
                                           shrinkage=c(.01,.001))                                           )
glmnet_learner = create.Learner("SL.glmnet",tune=list(alpha = seq(0,1,.2)))
ranger_learner = create.Learner("SL.ranger",params=list(min.node.size=30),
                                tune=list(num.trees=c(500,2500),mtry=c(2,3,4)))
knn_learner = create.Learner("SL.kernelKnn",tune=list(k=c(5,10,50)))
sl.lib <- c(ranger_learner$names,xgboost_learner$names,glmnet_learner$names,"SL.mean","SL.glm","SL.gam","SL.earth")

#tmle
W = subset(a,select=c(high_school,married,employed,white,age,BMI,BPS))
fit <- tmle(Y=as.numeric(a$outcome=="live birth"),
            A=a$treatment,
            W=data.frame(W),
            Q.SL.library=sl.lib, 
            g.SL.library=sl.lib,
            family="binomial")

g_pred <- fit$g$g1W

Q_pred <- fit$Qstar[,1]*(1-a$treatment) + fit$Qstar[,2]*(a$treatment)

Q1_pred <- fit$Qstar[,2] 

Q0_pred <- fit$Qstar[,1] 

mu0 <- Q0_pred
mu1 <- Q1_pred
mu <- Q_pred
pi <- g_pred

exposure <- a$treatment
outcome <- as.numeric(a$outcome=="live birth")
  
aipw <- function(mu0,mu1,mu,pi,exposure,outcome){
  ## risk difference
  aipw_RD <- mean((((2*exposure-1)*(outcome - mu))/((2*exposure-1)*pi + (1-exposure)) + mu1 - mu0))
  se_RD <- sd((((2*exposure-1)*(outcome - mu))/((2*exposure-1)*pi + (1-exposure)) + mu1 - mu0))/sqrt(length(exposure))
  aipw_RD.lcl <- aipw_RD - 1.96*se_RD
  aipw_RD.ucl <- aipw_RD + 1.96*se_RD
  
  aipw_eif1 <- (as.numeric(exposure==1)/pi)*(outcome - mu) + mu1
  aipw_eif0 <- (as.numeric(exposure==0)/pi)*(outcome - mu) + mu0
  
  sigma_covar <- matrix(data.frame(var(aipw_eif0),cov(aipw_eif0,aipw_eif1),cov(aipw_eif1,aipw_eif0),var(aipw_eif1)),nrow=2)
  
  ## risk ratio
  aipw_RR <- mean(aipw_eif1)/mean(aipw_eif0)
  se_RR <- ((as.numeric(sigma_covar[1,1])/(mean(aipw_eif0)^2)) - 
              (2*as.numeric(sigma_covar[1,2])/(mean(aipw_eif1)*mean(aipw_eif0))) + 
              (as.numeric(sigma_covar[2,2])/mean(aipw_eif1)^2) - 
              (2*as.numeric(sigma_covar[1,2])/(mean(aipw_eif1)*mean(aipw_eif0))))/sqrt(length(exposure))
  aipw_RR.lcl <- exp(log(aipw_RR) - 1.96*se_RR)
  aipw_RR.ucl <- exp(log(aipw_RR) + 1.96*se_RR)
  
  ## odds ratio
  aipw_OR <- (mean(aipw_eif1)/(1-mean(aipw_eif1))) / (mean(aipw_eif0)/(1-mean(aipw_eif0)))
  se_OR <- ((as.numeric(sigma_covar[1,1])/((mean(aipw_eif0)^2)*(mean(1-aipw_eif0)^2))) - 
              (2*as.numeric(sigma_covar[1,2])/(mean(aipw_eif1)*mean(aipw_eif0)*mean(1-aipw_eif1)*mean(1-aipw_eif0))) + 
              (as.numeric(sigma_covar[2,2])/((mean(aipw_eif1)^2)*(mean(1-aipw_eif1)^2))) - 
              (2*as.numeric(sigma_covar[1,2])/(mean(aipw_eif1)*mean(aipw_eif0)*mean(1-aipw_eif1)*mean(1-aipw_eif0))))/sqrt(length(exposure))
  aipw_OR.lcl <- exp(log(aipw_OR) - 1.96*se_OR)
  aipw_OR.ucl <- exp(log(aipw_OR) + 1.96*se_OR)
  
  res <- matrix(data.frame(aipw_RD,aipw_RD.lcl,aipw_RD.ucl,aipw_RR,aipw_RR.lcl,aipw_RR.ucl,aipw_OR,aipw_OR.lcl,aipw_OR.ucl),nrow=3,byrow=T)
  colnames(res) <- c("Estimate","95% LCL","95% UCL")
  row.names(res) <- c("Risk Difference","Risk Ratio", "Odds Ratio")
  return(res)
}

aipw(mu0=mu0,mu1=mu1,mu=mu,pi=pi,exposure=exposure,outcome=outcome)









