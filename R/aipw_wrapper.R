#' @title AIPW wrapper function
#'
#' @description
#' A wrapper function for `AIPW$new()$fit()$summary()`
#'
#' @param Y Outcome (binary integer: 0 or 1)
#' @param A Exposure (binary integer: 0 or 1)
#' @param verbose Whether to print the result (logical; Default = FALSE)
#' @param W covariates for both exposure and outcome models  (vector, matrix or data.frame). If null, this function will seek for
#' inputs from `W.Q` and `W.g`.
#' @param W.Q Only valid when `W` is null, otherwise it would be replaced by `W`.
#' Covariates for outcome model (vector, matrix or data.frame).
#' @param W.g Only valid when `W` is null, otherwise it would be replaced by `W`.
#' Covariates for exposure model (vector, matrix or data.frame)
#' @param Q.SL.library SuperLearner libraries or sl3 learner object (Lrnr_base) for outcome model
#' @param g.SL.library SuperLearner libraries or sl3 learner object (Lrnr_base) for exposure model
#' @param k_split Number of splitting (integer; range: from 1 to number of observation-1):
#'   if k_split=1, no cross-fitting;
#'   if k_split>=2, cross-fitting is used
#'                         (e.g., `k_split=10`, use 9/10 of the data to estimate and the remaining 1/10 leftover to predict).
#'   NOTE: it's recommended to use cross-fitting.
#' @param g.bound Value between \[0,1\] at which the propensity score should be truncated. Defaults to 0.025.
#' @param stratified_fit An indicator for whether the outcome model is fitted stratified by exposure status in the`fit()` method.
#'    Only when using `stratified_fit()` to turn on `stratified_fit = TRUE`, `summary` outputs average treatment effects among the treated and the controls.
#'
#' @export
#' @seealso [AIPW]
#' @return A fitted `AIPW` object with summarised results
#'
#' @examples
#' library(SuperLearner)
#' aipw_sl <- aipw_wrapper(Y=rbinom(100,1,0.5), A=rbinom(100,1,0.5),
#'                     W.Q=rbinom(100,1,0.5), W.g=rbinom(100,1,0.5),
#'                     Q.SL.library="SL.mean",g.SL.library="SL.mean",
#'                     k_split=1,verbose=FALSE)
aipw_wrapper = function(Y, A, verbose=TRUE,
                 W=NULL, W.Q=NULL, W.g=NULL,
                 Q.SL.library, g.SL.library,
                 k_split=10, g.bound=0.025,stratified_fit=FALSE){
  aipw_obj <- AIPW$new(Y=Y,A=A,verbose=verbose,
                       W=W, W.Q=W.Q,W.g=W.g,
                       Q.SL.library=Q.SL.library, g.SL.library=g.SL.library,
                       k_split=k_split)
  if (stratified_fit){
    aipw_obj$stratified_fit()
  } else{
    aipw_obj$fit()
  }
  aipw_obj$summary(g.bound=g.bound)

  invisible(aipw_obj)
}
