#' @title Augmented Inverse Probablity Weighting (AIPW) take TMLE or tmle3 inputs
#'
#' @description Define an R6Class aipw_tmle object whcih takes TMLE or tmle3 inputs
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @details create an AIPW object
#'
#' @return \code{AIPW} object
#'
#' @format \code{\link{R6Class}} object.
AIPW_tmle <- R6::R6Class(
  "AIPW_tmle",
  portable = TRUE,
  inherit = AIPW_base,
  public = list(
    #' @description
    #' Create a new `AIPW_tmle` object.
    #'
    #' @param Y outcome (binary integer: 0 or 1)
    #' @param A exposure (binary integer: 0 or 1)
    #' @param verbose whether to show progression bar and print the result (logical; Default = FALSE)
    #' @param tmle_fit a fitted `tmle` object
    #'
    #' @return A new `AIPW_tmle` obejct
    initialize = function(Y=NULL,A=NULL,tmle_fit = NULL,verbose=TRUE){
      #initialize from AIPW_base class
      super$initialize(Y=Y,A=A,verbose=verbose)
      #check the fitted object is tmle or tmle3 and import values accordingly
      if (any(class(tmle_fit) %in% "tmle")){
        message("Sample splitting was not supported with a fitted tmle object")
        self$obs_est$mu0 <- tmle_fit$Qstar[,1]
        self$obs_est$mu1 <- tmle_fit$Qstar[,2]
        self$obs_est$mu <- self$obs_est$mu0*(1-private$A) + self$obs_est$mu1*(private$A)
        self$obs_est$raw_p_score <- tmle_fit$g$g1W
      } else if (any(class(tmle_fit) %in% "tmle3_Fit")){
        message("Propensity scores from fitted tmle3 object are by default truncated (0.025)")
        #Q model
        Q.cf_task0 <- tmle_fit$tmle_params[[1]]$cf_likelihood$cf_tasks[[1]]
        Q.cf_task1 <- tmle_fit$tmle_params[[2]]$cf_likelihood$cf_tasks[[1]]
        self$obs_est$mu0 <- as.numeric(tmle_fit$likelihood$get_likelihood(Q.cf_task0,"Y", "validation"))
        self$obs_est$mu1 <- as.numeric(tmle_fit$likelihood$get_likelihood(Q.cf_task1,"Y", "validation"))
        self$obs_est$mu <-  self$obs_est$mu0*(1-private$A) + self$obs_est$mu1*(private$A)
        #g model
        #this p_score has already been truncated with 0.025=g.bound
        g.cf_task1 <- tmle_fit$tmle_task$generate_counterfactual_task("cf1",data.frame(A=1))
        self$obs_est$raw_p_score <- tmle_fit$likelihood$get_likelihood(g.cf_task1,"A")
      } else {
        stop("The tmle_fit is neither a `tmle` or `tmle3_Fit` object")
      }
    }
  )
)
