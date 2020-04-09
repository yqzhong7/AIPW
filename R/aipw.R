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
  #outcome model
  mu0 <- tmle_fit$Qstar[,1] #Q0_pred
  mu1 <- tmle_fit$Qstar[,2] #Q1_pred
  mu <- tmle_fit$Qstar[,1]*(1-exposure) + tmle_fit$Qstar[,2]*(exposure) #Q_pred
  #exposure model
  pi <- tmle_fit$g$g1W #g_pred (propensity score)

  ## risk difference
  aipw_i_est <- ((2*exposure-1)*(outcome - mu))/((2*exposure-1)*pi + (1-exposure)) + mu1 - mu0
  Z_norm <- sqrt(length(exposure))
  aipw_RD <- mean(aipw_i_est)
  se_RD <- stats::sd(aipw_i_est)/Z_norm
  aipw_RD.lcl <- aipw_RD - 1.96*se_RD
  aipw_RD.ucl <- aipw_RD + 1.96*se_RD
  aipw_RD.ci <- ci(aipw_RD,se_RD,ratio=F)

  ## risk ratio
  aipw_eif1 <- (as.numeric(exposure==1)/pi)*(outcome - mu) + mu1
  aipw_eif0 <- (as.numeric(exposure==0)/pi)*(outcome - mu) + mu0

  sigma_covar <- matrix(c(stats::var(aipw_eif0),
                                   stats::cov(aipw_eif0,aipw_eif1),
                                   stats::cov(aipw_eif1,aipw_eif0),
                                   stats::var(aipw_eif1)),nrow=2)


  aipw_RR <- mean(aipw_eif1)/mean(aipw_eif0)
  se_RR <- ((sigma_covar[1,1]/(mean(aipw_eif0)^2)) -
              (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0))) +
              (sigma_covar[2,2]/mean(aipw_eif1)^2) -
              (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0))))/Z_norm
  aipw_RR.ci <- ci(aipw_RR,se_RR,ratio=T)

  ## odds ratio
  aipw_OR <- (mean(aipw_eif1)/(1-mean(aipw_eif1))) / (mean(aipw_eif0)/(1-mean(aipw_eif0)))
  se_OR <- ((sigma_covar[1,1]/((mean(aipw_eif0)^2)*(mean(1-aipw_eif0)^2))) -
              (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0)*mean(1-aipw_eif1)*mean(1-aipw_eif0))) +
              (sigma_covar[2,2]/((mean(aipw_eif1)^2)*(mean(1-aipw_eif1)^2))) -
              (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0)*mean(1-aipw_eif1)*mean(1-aipw_eif0))))/Z_norm
  aipw_OR.ci <- ci(aipw_OR,se_OR,ratio=T)

  res <- matrix(c(aipw_RD,aipw_RD.ci,aipw_RR,aipw_RR.ci,aipw_OR,aipw_OR.ci),nrow=3,byrow=T)
  colnames(res) <- c("Estimate","95% LCL","95% UCL")
  row.names(res) <- c("Risk Difference","Risk Ratio", "Odds Ratio")
  return(res)
}
