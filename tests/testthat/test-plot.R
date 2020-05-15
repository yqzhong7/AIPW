#' @title Tesing plot.p_score
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/15
test_that("plot.p_score", {
  require(SuperLearner)
  require(ggplot2)
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
  expect_error(aipw$plot.p_score(),regexp = "Propensity scores are not estimated.")
  #after fitted
  aipw$fit()
  expect_message(g.plot <- aipw$plot.p_score(),regexp = "ATE has not been calculated.")
  expect_true(inherits(g.plot, "ggplot"))
  #after truncation
  capture.output(aipw$calculate_result())
  expect_silent(g.plot <- aipw$plot.p_score())
  expect_true(inherits(g.plot, "ggplot"))
})

