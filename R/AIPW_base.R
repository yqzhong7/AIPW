#' @title Augmented Inverse Probability Weighting Base Class (AIPW_base)
#'
#' @description A base class for AIPW that implements the common methods, such as \code{summary()} and \code{plot.p_score()}, inheritted by [AIPW] and [AIPW_tmle] class
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @return \code{AIPW} base object
#' @seealso [AIPW] and [AIPW_tmle]
#' @format \code{\link{R6Class}} object.
#' @export
AIPW_base <- R6::R6Class(
  "AIPW_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    #-------------------------public fields-----------------------------#
    #Number of observations
    n = NULL,
    #Components for estimating the influence functions of all observations to calculate average causal effects
    obs_est = list(mu0 = NULL,
                   mu1 = NULL,
                   mu = NULL,
                   raw_p_score = NULL,
                   p_score = NULL,
                   aipw_eif1 = NULL,
                   aipw_eif0 = NULL),
    #Risk difference, risk ratio, odds ratio and variance-covariance matrix for SE calculation
    estimates = list(RD = NULL,
                     RR = NULL,
                     OR = NULL,
                     sigma_covar = NULL),
    #A matrix contains RD, RR and OR with their SE and 95%CI
    result = NULL,
    #A density plot of propensity scores by exposure status (`ggplot2::geom_density`)
    g.plot = NULL,

    #-------------------------constructor-----------------------------#
    initialize = function(Y=NULL, A=NULL,verbose=TRUE){
      #save input into private fields
      private$Y=Y
      private$A=A
      private$verbose=verbose
      #check data length
      if (!(length(private$Y)==length(private$A))){
        stop("Please check the dimension of the data")
      }
      #setup
      self$n <- length(private$A)
      self$obs_est$mu0 <- rep(NA,self$n)
      self$obs_est$mu1 <- rep(NA,self$n)
      self$obs_est$mu <- rep(NA,self$n)
      self$obs_est$raw_p_score <- rep(NA,self$n)
    },

    #-------------------------summary method-----------------------------#
    summary = function(g.bound=0.025){
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
      self$obs_est$aipw_eif0 <- (as.numeric(private$A==0)/(1-self$obs_est$p_score))*(private$Y - self$obs_est$mu) + self$obs_est$mu0

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

    #-------------------------plot.p_score method-----------------------------#
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
      self$g.plot =  ggplot2::ggplot(data = plot_data,ggplot2::aes(x = p_score, group = A, color = A, fill=A)) +
        ggplot2::geom_density(alpha=0.5) +
        ggplot2::scale_x_continuous(limits = c(0,1)) +
        ggplot2::facet_wrap(~trunc) +
        ggtitle("Propensity scores by exposure status") +
        theme_bw()
      print(self$g.plot)
      invisible(self)
    }
  ),

  #-------------------------private fields and methods----------------------------#
  private = list(
    #input
    Y=NULL,
    A=NULL,
    verbose=NULL,
    g.bound=NULL,
    #private methods
    #Use individual estimates of efficient influence functions (obs_est$aipw_eif0 & obs_est$aipw_eif0) to calculate RD, RR and OR with SE and 95CI%
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



#' @name summary
#' @aliases summary.AIPW_base
#' @title Summary of the average treatment effects from AIPW
#'
#' @description
#' Calculate average causal effects in RD, RR and OR in the fitted [AIPW] or [AIPW_tmle] object using the estimated efficient influence functions
#'
#' @section R6 Usage:
#' \code{$summary(g.bound = 0.025)}
#'
#' @param g.bound Value between \[0,1\] at which the propensity score should be truncated. Defaults to 0.025.
#'
#' @seealso [AIPW] and [AIPW_tmle]
#'
#' @return `estimates` and `result`(public variables): Average treatment effect in RD, RR and OR
#'
#' @examples
#' library(SuperLearner)
#' aipw_sl <- AIPW$new(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
#'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
#'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
#'                     k_split=1,verbose=FALSE)$fit()
#' aipw_sl$summary(g.bound=0.025)
NULL



#' @name plot.p_score
#' @title Plot the propensity scores by exposure status
#'
#' @description
#' Plot and check the balance of propensity scores by exposure status
#'
#' @section R6 Usage:
#' \code{$plot.p_plot()}
#'
#' @seealso [AIPW] and [AIPW_tmle]
#'
#' @return `g.plot` (public variable): A density plot of propensity scores by exposure status (`ggplot2::geom_density`)
#'
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
#' aipw_sl$summary(g.bound=0.025)$plot.p_score()
NULL
