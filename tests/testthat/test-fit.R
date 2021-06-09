#' @title Testing fit: SuperLeaner & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/25
test_that("AIPW fit: SuperLeaner & k_split", {
  require(SuperLearner)
  ##k_split == 1: no cross-fitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE,
                    save.sl.fit = TRUE)
  aipw$fit()
  #check any null values after calculating results
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))

  ##k_split >=3: cross-fitting == k_split
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 2,verbose = FALSE)


  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,verbose = FALSE,
                    save.sl.fit = TRUE)
  aipw$fit()
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))
})

#' @title Testing fit: verbose and progressr
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/01/20
test_that("AIPW fit: verbose", {
  #verbose == TRUE: "Done!"
  library(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = T)
  expect_message(aipw$fit(),regexp = "Done!")

  ##progressr
  #when progressr not loaded
  expect_false(aipw$.__enclos_env__$private$isLoaded_progressr)
  #when loaded
  library(progressr)
  aipw$fit()
  expect_true(aipw$.__enclos_env__$private$isLoaded_progressr)
})

#' @title Testing fit: missing outcome reporting (N)
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/24
test_that("AIPW fit: missing outcome", {
  require(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  ##k_split == 1: no cross-fitting
  expect_warning(aipw <-  AIPW$new(Y=c(NA,vec()[2:100]),
                                   A=c(1,vec()[2:100]),
                                   W =vec(),
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE)$fit())
  #Check nuisance functions with missing data
  expect_true(is.na(aipw$obs_est$mu0[1]))
  expect_true(is.na(aipw$obs_est$mu1[1]))
  expect_true(is.na(aipw$obs_est$mu[1]))
  expect_false(is.na(aipw$obs_est$raw_p_score[1]))

  ##k_split == 2: 2 fold cross-fitting
  expect_warning(aipw <-  AIPW$new(Y=c(NA,vec()[2:100]),
                                   A=c(1,vec()[2:100]),
                                   W =vec(),
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 2,verbose = FALSE)$fit())
  #Check nuisance functions with missing data
  expect_true(is.na(aipw$obs_est$mu0[1]))
  expect_true(is.na(aipw$obs_est$mu1[1]))
  expect_true(is.na(aipw$obs_est$mu[1]))
  expect_false(is.na(aipw$obs_est$raw_p_score[1]))

  ##k_split == 3: 3 fold cross-fitting
  expect_warning(aipw <-  AIPW$new(Y=c(NA,vec()[2:100]),
                                   A=c(1,vec()[2:100]),
                                   W =vec(),
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 3,verbose = FALSE)$fit())
  #Check nuisance functions with missing data
  expect_true(is.na(aipw$obs_est$mu0[1]))
  expect_true(is.na(aipw$obs_est$mu1[1]))
  expect_true(is.na(aipw$obs_est$mu[1]))
  expect_false(is.na(aipw$obs_est$raw_p_score[1]))
})
