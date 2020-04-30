#' Title Get 95% Condifence Intervals
#'
#' @param est point estimate
#' @param se standard error
#' @param ratio logical (default==F); when TRUE, return exp(log(ci))
#'
#' @return lower and upper bounds of the 95% confidence interval
#'
#' @noRd
get_ci <- function(est, se, ratio=F) {
  if (ratio){
    est <- log(est)
    lcl <- est - 1.96*se
    ucl <- est + 1.96*se
    output <- exp(c(lcl,ucl))
  } else{
    lcl <- est - 1.96*se
    ucl <- est + 1.96*se
    output <- c(lcl,ucl)
  }
  return(output)
}
