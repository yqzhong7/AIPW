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
    initialize = function(Y=NULL,A=NULL,
                          tmle_fit = NULL){
      #initialize from AIPW_base class
      super$initialize(Y=Y,A=A)
      #check the fitted object is tmle or tmle3 and import values accordingly
      if (class(tmle_fitted)=="tmle"){
        self$obs_est$mu0 <- tmle_fit$Qstar[,1]
        self$obs_est$mu1 <- tmle_fit$Qstar[,2]
        self$obs_est$mu <- tmle_fit$Qstar[,1]*(1-private$A) + tmle_fit$Qstar[,2]*(private$A)
        self$obs_est$raw_p_score <- tmle_fit$g$g1W
      } else if (class(tmle_fitted)=="tmle3_Fit"){
        NULL
      } else {
        stop("The tmle_fitted is neither a `tmle` or `tmle3_Fit` object")
      }
    }
  )
)
