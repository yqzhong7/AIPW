#' @title Testing summary: SuperLeaner & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/08/09
test_that("AIPW summary: SuperLeaner & k_split", {
  require(SuperLearner)
  ##k_split == 1: no sample splitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = TRUE,
                    save.sl.fit = TRUE)
  #correctly print output
  expect_output(aipw$fit()$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))

  ##k_split >0: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,verbose = TRUE,
                    save.sl.fit = TRUE)
  #correctly print output
  expect_output(aipw$fit()$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))
})

#' @title Testing summary: sl3 & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/20
test_that("AIPW summary: sl3 & k_split", {
  require(sl3)
  ##k_split == 1: no sample splitting
  vec <- function() sample(0:1,100,replace = T)
  lrnr_glm <- sl3::Lrnr_glm$new()
  lrnr_mean <- sl3::Lrnr_mean$new()
  stacklearner <- sl3::Stack$new(lrnr_glm, lrnr_mean)
  metalearner <- sl3::Lrnr_nnls$new()
  sl3.lib <- sl3::Lrnr_sl$new(learners = stacklearner,
                              metalearner = metalearner)
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl3.lib,
                    g.SL.library=sl3.lib,
                    k_split = 1,verbose = TRUE,
                    save.sl.fit = TRUE)
  #correctly print output
  expect_output(aipw$fit()$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))
  #check results with RR and OR
  expect_equal(nrow(aipw$result),7)

  ##k_split >0: sample splitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl3.lib,
                    g.SL.library=sl3.lib,
                    k_split = 3,verbose = TRUE,
                    save.sl.fit = TRUE)
  #correctly print output
  expect_output(aipw$fit()$summary(), regexp = "Estimate")
  #check any null values after calculating results
  expect_false(any(sapply(aipw$estimates, is.null)))
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est, is.null)))
  expect_false(is.null(aipw$result))
  #check results with RR and OR
  expect_equal(nrow(aipw$result),7)
})


#' @title Testing summary: verbose
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/22
test_that("AIPW summary: verbose", {
  require(SuperLearner)
  ##k_split == 1: no sample splitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  ##verbose==False
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  expect_equal(capture.output(aipw$fit()$summary()),character(0))
})

#' @title Testing summary: continuous outcome reporting
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/22
test_that("AIPW summary: continuous outcome", {
  require(SuperLearner)
  ##k_split == 1: no sample splitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  ##verbose==False
  aipw <-  AIPW$new(Y=rnorm(100,10),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)$fit()$summary()
  #check results with RR and OR
  expect_equal(nrow(aipw$result),5)
})
