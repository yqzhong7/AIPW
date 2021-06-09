#' @title Testing stratified_fit: SuperLeaner & k_split
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/10
test_that("AIPW stratified_fit: SuperLeaner & k_split", {
  require(SuperLearner)
  ##k_split == 1: no cross-fitting
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE,
                    save.sl.fit = TRUE)
  expect_warning(aipw$stratified_fit())
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
  expect_warning(aipw$stratified_fit())
  expect_false(any(sapply(aipw$libs, is.null)))
  expect_false(any(sapply(aipw$obs_est[1:4], is.na))) #mu - raw_p_score
  expect_true(any(sapply(aipw$obs_est[5:7], is.null)))
  expect_true(is.null(aipw$result))
  expect_true(is.null(aipw$estimate))
})



#' @title Testing stratified_fit: verbose and progressr
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/10
test_that("AIPW stratified_fit: verbose", {
  #verbose == TRUE: "Done!"
  library(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = T)
  expect_message(aipw$stratified_fit(),regexp = "Done!")

  ##progressr
  #when loaded
  library(progressr)
  expect_message(aipw$stratified_fit())
  expect_true(aipw$.__enclos_env__$private$isLoaded_progressr)
})

#' @title Testing stratified_fit: missing outcome reporting (N)
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/24
test_that("AIPW stratified_fit: missing outcome", {
  require(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean")
  ##k_split == 1: no cross-fitting
  expect_warning(aipw <-  AIPW$new(Y=c(NA,vec()[2:100]),
                                   A=c(1,vec()[2:100]),
                                   W =vec(),
                                   Q.SL.library=sl.lib,
                                   g.SL.library=sl.lib,
                                   k_split = 1,verbose = FALSE)$stratified_fit())
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
                                   k_split = 2,verbose = FALSE)$stratified_fit())
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
                                   k_split = 3,verbose = FALSE)$stratified_fit())
  #Check nuisance functions with missing data
  expect_true(is.na(aipw$obs_est$mu0[1]))
  expect_true(is.na(aipw$obs_est$mu1[1]))
  expect_true(is.na(aipw$obs_est$mu[1]))
  expect_false(is.na(aipw$obs_est$raw_p_score[1]))
})


#' @title Testing stratified_fit: object
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2021/03/10
test_that("AIPW stratified_fit: object", {
  library(SuperLearner)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,
                    verbose = F,
                    save.sl.fit = T)
  # not stratified
  aipw$fit()
  expect_false(aipw$stratified_fitted)
  expect_equal(length(aipw$libs$Q.fit),3)


  #stratified
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 3,
                    verbose = F,
                    save.sl.fit = T)

  aipw$stratified_fit()
  expect_true(aipw$stratified_fitted)
  expect_equal(length(aipw$libs$Q.fit),3)
  expect_equal(length(aipw$libs$Q.fit[[1]]),2)
})
