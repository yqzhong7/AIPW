#' AIPW inputs
#'
#' Creating a input matrix for [`aipw()`]Augmented Inverse Probability Weighting
#'
#' @param Y a vector of outcome
#' @param A a vector of treatment
#' @param W a matrix or data.frame of covariates
#' @param Q.SL.library outcome model
#' @param g.SL.library expousre model (propensity score)
#' @param k_split number of splitting: if k_split=1, no sample splitting;
#'   if k_split>1, use similar technique of cross-validation
#'   (e.g., k_split=5, use 4/5 of the data to estimate the 1/5 leftover data)
#' @param verbose logical. Default = FALSE; when TRUE, show a text progress bar.
#'
#' @return a matrix with estimates for each row when treated aipw_eif1 and untreated aipw_eif0
#' @export
#'
aipw_input <- function(Y,A,W,Q.SL.library,g.SL.library,k_split=1,verbose=FALSE){
  #check data length
  if (!(length(Y)==length(A) & length(Y)==nrow(W))){
    stop("Please check the dimension of the data")
  }
  #setup
  n <- length(Y)
  mu0 <- rep(NA,n)
  mu1 <- rep(NA,n)
  mu <- rep(NA,n)
  pi <- rep(NA,n)
  k_index <- sample(rep(1:k_split,ceiling(n/k_split))[1:n],replace = F)

  pb = utils::txtProgressBar(min = 0, max = length(k_split), initial = 0,style = 3)
  for (i in 1:k_split){

    if (verbose){
      utils::setTxtProgressBar(pb,i)
    }

    #sample splitting
    if (k_split==1){
      train_index <- validation_index <- k_index==i
    } else{
      train_index <- k_index!=i
      validation_index <- k_index==i
    }
    train_set <- data.frame(cbind(A,W))[train_index,]
    validation_set <- data.frame(cbind(A,W))[validation_index,]

    #Q model(outcome model: g-comp)
    #fit with train set
    Q_fit <- SuperLearner::SuperLearner(Y = Y[train_index],
                          X = train_set,
                          SL.library = Q.SL.library,
                          family="binomial")
    #predict on validation set
    mu0[validation_index] <- as.numeric(stats::predict(Q_fit,newdata = transform(validation_set,A=0))$pred) #Q0_pred
    mu1[validation_index]  <- as.numeric(stats::predict(Q_fit,newdata = transform(validation_set,A=1))$pred) #Q1_pred
    mu[validation_index]  <- mu0[validation_index]*(1-A[validation_index]) + mu1[validation_index]*(A[validation_index]) #Q_pred

    #g model(exposure model: propensity score)
    #fit with train set
    g_fit <- SuperLearner::SuperLearner(Y=A[train_index],
                          X=W[train_index,],
                          SL.library = g.SL.library,
                          family="binomial")
    #predict on validation set
    pi[validation_index]  <- as.numeric(stats::predict(g_fit,newdata = W[validation_index,])$pred)  #g_pred
  }

  #AIPW est
  aipw_eif1 <- (as.numeric(A==1)/pi)*(Y - mu) + mu1
  aipw_eif0 <- (as.numeric(A==0)/pi)*(Y - mu) + mu0

  aipw_input_value <- matrix(c(aipw_eif1,aipw_eif0),ncol=2)
  colnames(aipw_input_value) <- c("aipw_eif1","aipw_eif0")
  return(aipw_input_value)
}


#' AIPW
#'
#' Augmented Inverse Probability Weighting
#'
#' @param aipw_input a matrix of individual estimates created by  [`aipw_input()`]
#' @param A exposure (only needed when using tmle fitted object as input)
#' @param Y outcome (only needed when using tmle fitted object as input)
#' @param tmle_fit fitted `tmle`` object
#'
#' @return a matrix with risk difference (RD), risk ratio (RR), odds ratio (OR)
#' @export
#'
aipw <- function(aipw_input=NULL,tmle_fit=NULL,A=NULL,Y=NULL){
  if (is.null(aipw_input) & !is.null(tmle_fit)){
    if (class(tmle_fit)!="tmle") {
      stop("Input for 'tmle_fit' is an invalid TMLE fitted object")
    }
    #extract from tmle object and convert to aipw input
    #outcome model
    mu0 <- tmle_fit$Qstar[,1] #Q0_pred
    mu1 <- tmle_fit$Qstar[,2] #Q1_pred
    mu <- tmle_fit$Qstar[,1]*(1-A) + tmle_fit$Qstar[,2]*(A) #Q_pred
    #exposure model
    pi <- tmle_fit$g$g1W #g_pred (propensity score)
    #AIPW est
    aipw_eif1 <- (as.numeric(A==1)/pi)*(Y - mu) + mu1
    aipw_eif0 <- (as.numeric(A==0)/pi)*(Y - mu) + mu0

  }

  if (!is.null(aipw_input) & is.null(tmle_fit)){
    aipw_eif1 <- aipw_input[,1]
    aipw_eif0 <- aipw_input[,2]
  }

  Z_norm <- sqrt(length(aipw_eif1))

  ## risk difference
  aipw_RD <- mean(aipw_eif1 - aipw_eif0)
  se_RD <- stats::sd(aipw_eif1 - aipw_eif0)/Z_norm
  aipw_RD.ci <- ci(aipw_RD,se_RD,ratio=F)

  ## risk ratio
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
              (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0)
                                   *mean(1-aipw_eif1)*mean(1-aipw_eif0))))/Z_norm
  aipw_OR.ci <- ci(aipw_OR,se_OR,ratio=T)

  N <- length(aipw_eif1)

  res <- matrix(c(aipw_RD,se_RD,aipw_RD.ci,N,aipw_RR,se_RR,aipw_RR.ci,N,aipw_OR,se_OR,aipw_OR.ci,N),nrow=3,byrow=T)
  colnames(res) <- c("Estimate","SE","95% LCL","95% UCL","N")
  row.names(res) <- c("Risk Difference","Risk Ratio", "Odds Ratio")
  return(res)
}
