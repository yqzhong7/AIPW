#' @title Augmented Inverse Probablity Weighting Base Class (AIPW_base)
#'
#' @description Define a base R6Class for later subclasses inheritance of
#' \code{AIPW_base$calculate_result()} and \code{AIPW_base$plot.p_score()}
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @return \code{AIPW} base object
#'
#' @format \code{\link{R6Class}} object.
AIPW_base <- R6::R6Class(
  "AIPW_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    #' @field n number of observations
    n = NULL,
    #' @field obs_est components for estimating the influence functions of all observations to calculate average causal effects
    obs_est = list(mu0 = NULL,
                   mu1 = NULL,
                   mu = NULL,
                   raw_p_score = NULL,
                   p_score = NULL,
                   aipw_eif1 = NULL,
                   aipw_eif0 = NULL),
    #' @field estimates risk difference, risk ratio, odds ratio and variance-covariance matrix for SE calculation
    estimates = list(RD = NULL,
                     RR = NULL,
                     OR = NULL,
                     sigma_covar = NULL),
    #' @field result a matrix contains RD, RR and OR with their SE and 95%CI
    result = NULL,

    #' @description
    #' Create a new `AIPW_base` object.
    #'
    #' @param Y outcome (binary integer: 0 or 1)
    #' @param A exposure (binary integer: 0 or 1)
    #' @param ... not used
    #' @return A new `AIPW_base` obejct
    #'
    #' @examples
    #' library(SuperLearner)
    #' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
    #'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
    #'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
    #'                     k_split=1,verbose=FALSE)
    initialize = function(Y=NULL, A=NULL){
      #save input into private fields
      private$Y=Y
      private$A=A
      #check data length
      if (!(length(private$Y)==length(private$A))){
        stop("Please check the dimension of the data")
      }
      #setup
      self$n <- length(private$Y)
      self$obs_est$mu0 <- rep(NA,self$n)
      self$obs_est$mu1 <- rep(NA,self$n)
      self$obs_est$mu <- rep(NA,self$n)
      self$obs_est$raw_p_score <- rep(NA,self$n)
    },
    #' @description
    #' Calculate average causal effects in RD, RR and OR in the fitted `AIPW` obejct with the estimated influence functions
    #'
    #' @param g.bound value between \[0,1\] at which the propensity score should be truncated. Defaults to 0.025.
    #'
    #' @return An `AIPW` obejct with average treatment effect estimations in RD, RR and OR
    #'
    #' @examples
    #' library(SuperLearner)
    #' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
    #'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
    #'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
    #'                     k_split=1,verbose=FALSE)$fit()
    #' aipw_sl$calculate_result(g.bound=0.025)
    calculate_result = function(g.bound=0.025){
      #p_score truncation
      private$g.bound=g.bound
      #check g.bound value
      if (!is.numeric(private$g.bound)){
        stop("g.bound must be a numeric value")
      } else if (private$g.bound>1|private$g.bound<0){
        stop("g.bound must between 0 and 1")
      }
      self$obs_est$p_score <- private$.bound(self$obs_est$raw_p_score)

      #AIPW est
      self$obs_est$aipw_eif1 <- (as.numeric(private$A==1)/self$obs_est$p_score)*(private$Y - self$obs_est$mu) + self$obs_est$mu1
      self$obs_est$aipw_eif0 <- (as.numeric(private$A==0)/self$obs_est$p_score)*(private$Y - self$obs_est$mu) + self$obs_est$mu0

      root_n <- sqrt(self$n)

      ## risk difference
      self$estimates$RD <- private$get_RD(self$obs_est$aipw_eif1, self$obs_est$aipw_eif0, root_n)

      ## var-cov mat for rr and or calculation
      self$estimates$sigma_covar <- private$get_sigma_covar(self$obs_est$aipw_eif0,self$obs_est$aipw_eif1)

      ## risk ratio
      self$estimates$RR <- private$get_RR(self$obs_est$aipw_eif1,self$obs_est$aipw_eif0, self$estimates$sigma_covar, root_n)

      ## odds ratio
      self$estimates$OR <- private$get_OR(self$obs_est$aipw_eif1,self$obs_est$aipw_eif0, self$estimates$sigma_covar, root_n)

      self$result <- cbind(matrix(c(self$estimates$RD,self$estimates$RR,self$estimates$OR),nrow=3,byrow=T),self$n)
      colnames(self$result) <- c("Estimate","SE","95% LCL","95% UCL","N")
      row.names(self$result) <- c("Risk Difference","Risk Ratio", "Odds Ratio")
      if (private$verbose){
        print(self$result,digit=3)
      }
      invisible(self)
    },
    #' @description
    #' Plot and check the balance of propensity scores by exposure status
    #'
    #' @return A density plot of propensity scores by exposure status (`ggplot2::geom_density`)
    #' @examples
    #' library(SuperLearner)
    #' library(ggplot2)
    #' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
    #'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
    #'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
    #'                     k_split=1,verbose=FALSE)$fit()
    #' #before average treatment effect calculation
    #' aipw_sl$plot.p_score()
    #' #after calculation
    #' aipw_sl$calculate_result(g.bound=0.025)$plot.p_score()
    plot.p_score = function(){
      #check if ggplot2 library is loaded
      if (!any(names(sessionInfo()$otherPkgs) %in% c("ggplot2"))){
        stop("`ggplot2` package is not loaded.")
      }
      #input check
      if (any(is.na(self$obs_est$raw_p_score))){
        stop("Propensity scores are not estimated.")
      } else if (is.null(self$obs_est$p_score)) {
        #p_score before truncation (estimated ps)
        plot_data = data.frame(A = factor(private$A),
                               p_score= self$obs_est$raw_p_score,
                               trunc = "Not truncated")
        message("ATE has not been calculated.")
      } else {
        plot_data = rbind(data.frame(A = factor(private$A),
                                     p_score= self$obs_est$raw_p_score,
                                     trunc = "Not truncated"),
                          data.frame(A = factor(private$A),
                                     p_score= self$obs_est$p_score,
                                     trunc = "Truncated"))
      }
      g.plot <-  ggplot2::ggplot(data = plot_data,ggplot2::aes(x = p_score, group = A, color = A, fill=A)) +
        ggplot2::geom_density(alpha=0.5) +
        ggplot2::scale_x_continuous(limits = c(0,1)) +
        ggplot2::facet_wrap(~trunc) +
        ggtitle("Propensity scores by exposure status") +
        theme_bw()
      print(g.plot)
    }
  ),
  private = list(
    #input
    Y=NULL,
    A=NULL,
    g.bound=NULL,
    #private methods
    #Use individaul estimates (obs_est$aipw_eif0 & obs_est$aipw_eif0 ) to calcualte RD, RR and OR with SE and 95CI%
    get_RD = function(aipw_eif1,aipw_eif0,root_n){
      est <- mean(aipw_eif1 - aipw_eif0)
      se <- stats::sd(aipw_eif1 - aipw_eif0)/root_n
      ci <- get_ci(est,se,ratio=F)
      output = c(est, se, ci)
      names(output) = c("Estimate","SE","95% LCL","95% UCL")
      return(output)
    },
    get_RR = function(aipw_eif1,aipw_eif0,sigma_covar,root_n){
      est <- mean(aipw_eif1)/mean(aipw_eif0)
      se <- sqrt((sigma_covar[1,1]/(mean(aipw_eif0)^2)) -
                   (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0))) +
                   (sigma_covar[2,2]/mean(aipw_eif1)^2) -
                   (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0))))/root_n
      ci <- get_ci(est,se,ratio=T)
      output = c(est, se, ci)
      names(output) = c("Estimate","SE","95% LCL","95% UCL")
      return(output)
    },
    get_OR = function(aipw_eif1,aipw_eif0,sigma_covar,root_n){
      est <- (mean(aipw_eif1)/(1-mean(aipw_eif1))) / (mean(aipw_eif0)/(1-mean(aipw_eif0)))
      se <- sqrt((sigma_covar[1,1]/((mean(aipw_eif0)^2)*(mean(1-aipw_eif0)^2))) -
                   (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0)*mean(1-aipw_eif1)*mean(1-aipw_eif0))) +
                   (sigma_covar[2,2]/((mean(aipw_eif1)^2)*(mean(1-aipw_eif1)^2))) -
                   (2*sigma_covar[1,2]/(mean(aipw_eif1)*mean(aipw_eif0)
                                        *mean(1-aipw_eif1)*mean(1-aipw_eif0))))/root_n
      ci <- get_ci(est,se,ratio=T)
      output = c(est, se, ci)
      names(output) = c("Estimate","SE","95% LCL","95% UCL")
      return(output)
    },
    get_sigma_covar = function(aipw_eif0,aipw_eif1){
      mat <- matrix(c(stats::var(aipw_eif0),
                      stats::cov(aipw_eif0,aipw_eif1),
                      stats::cov(aipw_eif1,aipw_eif0),
                      stats::var(aipw_eif1)),nrow=2)
      return(mat)
    },
    #setup the bounds for the propensity score to ensure the balance
    .bound = function(p_score,bound = private$g.bound){
      res <- base::ifelse(p_score<bound,bound,
                          base::ifelse(p_score>(1-bound),(1-bound),p_score))
      return(res)
    }
  )
)
