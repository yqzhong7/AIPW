#' AIPW
#'
#' Augmented Inverse Probability Weighting
#'
#' @param exposure binary numeric vector
#' @param outcome binary numeric vector
#' @param tmle_fit fitted tmle object
#' @param tmle3_fit place holder for tmle3 function
#'
#' @return a matrix with risk difference (RD), risk ratio (RR), odds ratio (OR)
#' @export
#'
aipw <- function(exposure,outcome,tmle_fit,tmle3_fit=NULL){
  if (class(tmle_fit)!="tmle") {
    stop("Input for 'tmle_fit' is an invalid TMLE fitted object")
    }

  #extract from tmle object and convert to aipw input
  mu0 <- tmle_fit$Qstar[,1]
  mu1 <- tmle_fit$Qstar[,2]
  mu <- tmle_fit$Qstar[,1]*(1-exposure) + tmle_fit$Qstar[,2]*(exposure)
  pi <- tmle_fit$g$g1W

  ## risk difference
  aipw_RD <- mean((((2*exposure-1)*(outcome - mu))/((2*exposure-1)*pi + (1-exposure)) + mu1 - mu0))
  se_RD <- stats::sd((((2*exposure-1)*(outcome - mu))/((2*exposure-1)*pi + (1-exposure)) + mu1 - mu0))/sqrt(length(exposure))
  aipw_RD.lcl <- aipw_RD - 1.96*se_RD
  aipw_RD.ucl <- aipw_RD + 1.96*se_RD

  aipw_eif1 <- (as.numeric(exposure==1)/pi)*(outcome - mu) + mu1
  aipw_eif0 <- (as.numeric(exposure==0)/pi)*(outcome - mu) + mu0

  sigma_covar <- matrix(data.frame(stats::var(aipw_eif0),
                                   stats::cov(aipw_eif0,aipw_eif1),
                                   stats::cov(aipw_eif1,aipw_eif0),
                                   stats::var(aipw_eif1)),nrow=2)

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
