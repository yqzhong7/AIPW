#' @title Testing AIPW_nuis: user input nuisance function
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/28
test_that("AIPW_nuis class", {
  vec <- function() sample(0:1,100,replace = T)
  con_vec <- function() {
    x <- rnorm(100, mean = 0.5, sd = 0.5/3)
    x <- pmax(0, x)
    x <- pmin(1, x)
    return(x)
  }
  #test constructor
  expect_message(aipw_nuis <- AIPW_nuis$new(A=vec(),Y=vec(),
                        mu0 = con_vec(), mu1 = con_vec(), raw_p_score = con_vec()),
                 info = "Cross-fitting for estimating nuisance functions is recommended")
  #correctly print output
  expect_output(aipw_nuis$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw_nuis$estimates, is.null)))
  expect_false(any(sapply(aipw_nuis$libs, is.null)))
  expect_false(any(sapply(aipw_nuis$obs_est, is.null)))
  expect_false(is.null(aipw_nuis$result))
})
