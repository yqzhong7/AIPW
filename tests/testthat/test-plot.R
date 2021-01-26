#' @title Testing plot.p_score
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/27
test_that("plot.p_score", {
  require(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  expect_error(aipw$plot.p_score(),regexp = "`ggplot2` package is not loaded.")
  require(ggplot2)
  #before fitted
  expect_error(aipw$plot.p_score(),regexp = "Propensity scores are not estimated.")
  #after fitted
  aipw$fit()
  expect_message(g.plot <- aipw$plot.p_score(),regexp = "ATE has not been calculated.")
  expect_true(inherits(aipw$g.plot, "ggplot"))
  #after truncation
  capture.output(aipw$fit()$summary())
  expect_silent(g.plot <- aipw$plot.p_score())
  expect_true(inherits(aipw$g.plot, "ggplot"))
})

#' @title Testing plot.ip_weights
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/26
test_that("plot.ip_weights", {
  require(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  #before fitted
  expect_error(aipw$plot.ip_weights(),regexp = "Propensity scores are not estimated.")
  #after fitted
  aipw$fit()
  expect_error(ip_weights.plot <- aipw$plot.ip_weights(),regexp = "ATE has not been calculated.")
  #after truncation
  capture.output(aipw$fit()$summary())
  expect_silent(ip_weights.plot <- aipw$plot.ip_weights())
  expect_true(inherits(aipw$ip_weights.plot, "ggplot"))
})

