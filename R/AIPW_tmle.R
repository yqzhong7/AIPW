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
    #' @param Y Outcome (binary integer: 0 or 1)
    #' @param A Exposure (binary integer: 0 or 1)
    #' @param verbose Whether to show progression bar and print the result (logical; Default = FALSE)
    #' @param tmle_fit A fitted `tmle` object
    #'
    #' @return A new `AIPW_tmle` obejct
    #'
    #' @examples
    #' \dontrun{
    #' vec <- function() sample(0:1,100,replace = TRUE)
    #' df <- data.frame(replicate(4,vec()))
    #' names(df) <- c("A","Y","W1","W2")
    #'
    #' ## From tmle
    #' require(tmle)
    #' require(SuperLearner)
    #' tmle_fit <- tmle(Y=df$Y,A=df$A,W=subset(df,select=c("W1","W2")),
    #'                  Q.SL.library="SL.glm",
    #'                  g.SL.library="SL.glm",
    #'                  family="binomial")
    #' AIPW_tmle$new(A=df$A,Y=df$Y,tmle_fit = tmle_fit,verbose = TRUE)$calculate_result()
    #'
    #'
    #' ## From tmle3
    #' # tmle3 simple implementation
    #' require(tmle3)
    #' require(sl3)
    #' node_list <- list(A = "A",Y = "Y",W = c("W1","W2"))
    #' or_spec <- tmle_OR(baseline_level = "0",contrast_level = "1")
    #' tmle_task <- or_spec$make_tmle_task(df,node_list)
    #' lrnr_glm <- make_learner(Lrnr_glm)
    #' sl <- Lrnr_sl$new(learners = list(lrnr_glm))
    #' learner_list <- list(A = sl, Y = sl)
    #' tmle3_fit <- tmle3(or_spec, data=df, node_list, learner_list)
    #'
    #' # parse tmle3_fit into AIPW_tmle class
    #' AIPW_tmle$new(A=df$A,Y=df$Y,tmle_fit = tmle3_fit,verbose = TRUE)$calculate_result()
    #' }
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
