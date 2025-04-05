#' @title Augmented Inverse Probability Weighting Base Class (AIPW_base)
#'
#' @description A base class for AIPW that implements the common methods, such as \code{summary()} and \code{plot.p_score()}, inheritted by [AIPW] and [AIPW_tmle] class
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom stats predict
#' @importFrom utils head
#' @importFrom ggplot2 ggplot aes geom_density theme_bw labs
#'
#' @return \code{AIPW} base object
#' @seealso [AIPW] and [AIPW_tmle]
#' @format \CRANpkg{R6} object.
#' @export
AIPW_base <- R6::R6Class(
  "AIPW_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    #-------------------------public fields-----------------------------#
    #Number of observations
    n = NULL,
    #Number of exposed
    n_A1 = NULL,
    #Number of unexposed
    n_A0 = NULL,
    #Fit the outcome model stratified by exposure status (only applicable to AIPW class or manual setup)
    stratified_fitted = FALSE,
    #Components for estimating the influence functions of all observations to calculate average causal effects
    obs_est = list(mu0 = NULL,
                   mu1 = NULL,
                   mu = NULL,
                   raw_p_score = NULL,
                   p_score = NULL,
                   ip_weights = NULL,
                   aipw_eif1 = NULL,
                   aipw_eif0 = NULL),
    #ATE: Risk difference, risk ratio, odds ratio and variance-covariance matrix for SE calculation
    estimates = list(risk_A1 = NULL,
                     risk_A0 = NULL,
                     RD = NULL,
                     RR = NULL,
                     OR = NULL,
                     sigma_covar = NULL),
    #ATT: Risk difference
    ATT_estimates = list(RD = NULL),
    #ATC: Risk difference
    ATC_estimates = list(RD = NULL),
    #A matrix contains RD, RR and OR with their SE and 95%CI
    result = NULL,
    #A density plot of propensity scores by exposure status (`ggplot2::geom_density`)
    g.plot = NULL,
    #A box plot of inverse probability weights using truncated propensity scores by exposure status (`ggplot2::geom_boxplot`)
    ip_weights.plot = NULL,

    #-------------------------constructor-----------------------------#
    initialize = function(Y=NULL, A=NULL,verbose=TRUE){
      #save input into private fields
      private$Y=as.numeric(Y)
      private$A=as.numeric(A)
      private$observed = as.numeric(!is.na(private$Y))
      private$verbose=verbose
      #check data length
      if (length(private$Y)!=length(private$A)){
        stop("Please check the dimension of the data")
      }
      #detect outcome is binary or continuous
      if (length(unique(private$Y[!is.na(private$Y)]))==2) {
        private$Y.type = 'binomial'
      } else {
        private$Y.type = 'gaussian'
      }

      #check missing exposure
      if (any(is.na(private$A))){
        stop("Missing exposure is not allowed.")
      }

      #check missing outcome
      if (any(private$observed == 0)){
        warning("Missing outcome is detected. Analysis assumes missing at random (MAR).")
        private$Y.missing = TRUE
      }

      #setup
      private$AxObserved = private$A * private$observed #I(A=a, observed==1)
      self$n <- length(private$A)
      self$n_A1 <- sum(private$A==1)
      self$n_A0 <- sum(private$A==0)
      self$obs_est$mu0 <- rep(NA,self$n)
      self$obs_est$mu1 <- rep(NA,self$n)
      self$obs_est$mu <- rep(NA,self$n)
      self$obs_est$raw_p_score <- rep(NA,self$n)
    },

    #-------------------------summary method-----------------------------#
    summary = function(g.bound=0.025){
      #p_score truncation
      if (length(g.bound) > 2){
        warning('More than two `g.bound` are provided. Only the first two will be used.')
        g.bound = g.bound[1:2]
      } else if (length(g.bound) ==1 & g.bound[1] >= 0.5){
          stop("`g.bound` >= 0.5 is not allowed when only one `g.bound` value is provided")
      }

      private$g.bound=g.bound
      #check g.bound value
      if (!is.numeric(private$g.bound)){
        stop("`g.bound` must be numeric")
      } else if (max(private$g.bound) > 1 | min(private$g.bound) < 0){
        stop("`g.bound` must between 0 and 1")
      }
      self$obs_est$p_score <- private$.bound(self$obs_est$raw_p_score)

      #inverse probability weights
      self$obs_est$ip_weights <- (as.numeric(private$A==1)/self$obs_est$p_score) + (as.numeric(private$A==0)/(1-self$obs_est$p_score))

      ##------AIPW est------##
      #### ATE EIF
      self$obs_est$aipw_eif1 <- ifelse(private$observed == 1,
                                       (as.numeric(private$A[private$observed==1]==1)/self$obs_est$p_score[private$observed==1])*
                                         (private$Y[private$observed==1] - self$obs_est$mu[private$observed==1]) +
                                         self$obs_est$mu1[private$observed==1],
                                       0)
      self$obs_est$aipw_eif0 <- ifelse(private$observed == 1,
                                       (as.numeric(private$A[private$observed==1]==0)/(1-self$obs_est$p_score[private$observed==1]))*
                                         (private$Y[private$observed==1] - self$obs_est$mu[private$observed==1]) +
                                         self$obs_est$mu0[private$observed==1],
                                       0)

      root_n <- sqrt(self$n)

      ## risk for the treated and controls
      self$estimates$risk_A1 <- private$get_RD(self$obs_est$aipw_eif1, 0, root_n)
      self$estimates$risk_A0 <- private$get_RD(self$obs_est$aipw_eif0, 0, root_n)

      ## risk difference
      self$estimates$RD <- private$get_RD(self$obs_est$aipw_eif1, self$obs_est$aipw_eif0, root_n)

      #results on additive scales
      self$result <- cbind(matrix(c(self$estimates$risk_A1, self$estimates$risk_A0,
                                    self$estimates$RD), nrow=3, byrow=T),
                           c( self$n_A1, self$n_A0,rep(self$n,1)))
      row.names(self$result) <- c("Risk of Exposure", "Risk of Control","Risk Difference")
      colnames(self$result) <- c("Estimate","SE","95% LCL","95% UCL","N")

      if (private$Y.type == 'binomial'){
        ## var-cov mat for rr and or calculation
        self$estimates$sigma_covar <- private$get_sigma_covar(self$obs_est$aipw_eif0,self$obs_est$aipw_eif1)

        ## risk ratio
        self$estimates$RR <- private$get_RR(self$obs_est$aipw_eif1,self$obs_est$aipw_eif0, self$estimates$sigma_covar, root_n)

        ## odds ratio
        self$estimates$OR <- private$get_OR(self$obs_est$aipw_eif1,self$obs_est$aipw_eif0, self$estimates$sigma_covar, root_n)

        #w/ results on the multiplicative scale
        mult_result <- cbind(matrix(c(self$estimates$RR, self$estimates$OR),nrow=2,byrow=T),self$n)
        row.names(mult_result) <- c("Risk Ratio", "Odds Ratio")
        self$result <- rbind(self$result, mult_result)
      }

      #### ATT/ATC
      if (self$stratified_fitted) {
        #ATT
        self$ATT_estimates$RD <- private$get_ATT_RD(mu0 = self$obs_est$mu0[private$observed==1],
                                       p_score = self$obs_est$p_score[private$observed==1],
                                       A_level = 1, root_n=root_n, ATC = F)
        self$ATC_estimates$RD <- private$get_ATT_RD(mu0 = self$obs_est$mu1[private$observed==1],
                                                    p_score = 1-self$obs_est$p_score[private$observed==1],
                                                    A_level = 0, root_n=root_n, ATC = T)
        ATT_ATC_result <- matrix(c(self$ATT_estimates$RD, self$n,
                                   self$ATC_estimates$RD, self$n), nrow = 2,byrow = T)
        row.names(ATT_ATC_result) <- c("ATT Risk Difference","ATC Risk Difference")
        self$result <- rbind(self$result, ATT_ATC_result)
      }

      #### Change row names for continuous outcome
      if (private$Y.type == 'gaussian'){
        row.names(self$result) = gsub("Risk", "Mean", row.names(self$result))
      }

      if (private$verbose){
        print(self$result,digit=3)
      }
      invisible(self)
    },

    #-------------------------plot.p_score method-----------------------------#
    plot.p_score = function(print.ip_weights = F){
      #check if ggplot2 library is loaded
      if (!any(names(sessionInfo()$otherPkgs) %in% c("ggplot2"))){
        stop("`ggplot2` package is not loaded.")
      }

      plot_data_A = factor(private$A, levels = 0:1)

      #input check
      if (any(is.na(self$obs_est$raw_p_score))){
        stop("Propensity scores are not estimated.")
      } else if (is.null(self$obs_est$p_score)) {
        #p_score before truncation (estimated ps)
        plot_data = data.frame(A = plot_data_A,
                               p_score= self$obs_est$raw_p_score,
                               trunc = "Not truncated")
        message("ATE has not been calculated.")
      } else {
        plot_data = rbind(data.frame(A = plot_data_A,
                                     p_score= self$obs_est$raw_p_score,
                                     trunc = "Not truncated"),
                          data.frame(A = plot_data_A,
                                     p_score= self$obs_est$p_score,
                                     trunc = "Truncated"))
      }
      self$g.plot =  ggplot2::ggplot(data = plot_data,ggplot2::aes(x = p_score, group = A, color = A, fill=A)) +
        ggplot2::geom_density(alpha=0.5) +
        ggplot2::scale_x_continuous(limits = c(0,1)) +
        ggplot2::facet_wrap(~trunc) +
        ggtitle("Propensity scores by exposure status") +
        theme_bw() +
        theme(legend.position = 'bottom')
        xlab('Propensity Scores')

      print(self$g.plot)
      invisible(self)
    }
  ,
  #-------------------------plot.ip_weights method-----------------------------#
  plot.ip_weights = function(){
    #check if ggplot2 library is loaded
    if (!any(names(sessionInfo()$otherPkgs) %in% c("ggplot2"))){
      stop("`ggplot2` package is not loaded.")
    }

    plot_data_A = factor(private$A, levels = 0:1)

    #input check
    if (any(is.na(self$obs_est$raw_p_score))){
      stop("Propensity scores are not estimated.")
    } else if (is.null(self$obs_est$p_score)) {
      stop("ATE has not been calculated.")
    } else {
      ipw_plot_data = data.frame(A = plot_data_A, ip_weights= self$obs_est$ip_weights)
      self$ip_weights.plot = ggplot2::ggplot(data = ipw_plot_data, ggplot2::aes(y = ip_weights, x = A, fill = A)) +
        ggplot2::geom_boxplot(alpha=0.5) +
        ggtitle("IP-weights using truncated propensity scores by exposure status") +
        theme_bw() +
        ylab('Inverse Probablity Weights') +
        coord_flip() +
        theme(legend.position = 'bottom')

      print(self$ip_weights.plot)
    }

    invisible(self)
  }
),
  #-------------------------private fields and methods----------------------------#
  private = list(
    #input
    Y=NULL,
    A=NULL,
    observed=NULL,
    AxObserved = NULL,
    verbose=NULL,
    g.bound=NULL,
    #outcome type
    Y.type = NULL,
    Y.missing = FALSE,
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
    #ATT/ATC calculation
    get_ATT_RD = function(A =private$A[private$observed==1], Y = private$Y[private$observed==1],
                          mu0, p_score, A_level, root_n, ATC = F){
      I_A = (A==A_level) / mean(A==A_level)
      I_A_com = (1-A==A_level) / mean(1-(A==A_level))
      eif <- I_A*Y  - (I_A*(mu0) + I_A_com*(Y-mu0)*p_score/(1-p_score))
      est <- mean(eif)
      if (ATC){
        est <- -1 * est
      }
      se <- stats::sd(eif - I_A*est)/root_n
      ci <- get_ci(est,se,ratio=F)
      output = c(est, se, ci)
      names(output) = c("Estimate","SE","95% LCL","95% UCL")
      return(output)
    },
    #setup the bounds for the propensity score to ensure the balance
    .bound = function(p_score,bound = private$g.bound){
      if (length(bound) == 1){
        res <- base::ifelse(p_score<bound, bound,
                            base::ifelse(p_score > (1-bound), (1-bound) ,p_score))
      } else {
        res <- base::ifelse(p_score< min(bound), min(bound),
                            base::ifelse(p_score > max(bound), max(bound), p_score))
      }
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
#' \code{$summary(g.bound = 0.025)} \cr
#' \code{$summary(g.bound = c(0.025,0.975))}
#'
#' @param g.bound Value between \[0,1\] at which the propensity score should be truncated.
#' Propensity score will be truncated to \eqn{[g.bound, 1-g.bound]} when one g.bound value is provided, or to \eqn{[min(g.bound), max(g.bound)]} when two values are provided.
#'  \strong{Defaults to 0.025}.
#'
#' @seealso [AIPW] and [AIPW_tmle]
#'
#' @return `estimates` and `result` (public variables): Risks, Average treatment effect in RD, RR and OR.
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
NULL

#' @name plot.ip_weights
#' @title Plot the inverse probability weights using truncated propensity scores by exposure status
#'
#' @description
#' Plot and check the balance of propensity scores by exposure status
#'
#' @section R6 Usage:
#' \code{$plot.ip_weights()}
#'
#' @seealso [AIPW] and [AIPW_tmle]
#'
#' @return `ip_weights.plot` (public variable): A box plot of inverse probability weights using truncated propensity scores by exposure status (`ggplot2::geom_boxplot`)
NULL
