#' @title Testing future_lapply wrapper
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/26
test_that("AIPW .flappy: lapply", {
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
  #check function body is lapply rather than future_lapply
  func_body<- as.character(body(aipw$.__enclos_env__$private$.f_lapply))
  expect_true(any(grepl("^lapply",func_body)))
})

test_that("AIPW .flappy: future_lapply", {
  require(SuperLearner)
  require(future.apply)
  vec <- function() sample(0:1,100,replace = T)
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec(),
                    A=vec(),
                    W.Q =vec(),
                    W.g =vec(),
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  #check function body is lapply rather than future_lapply
  func_body<- as.character(body(aipw$.__enclos_env__$private$.f_lapply))
  expect_true(any(grepl("future_lapply",func_body)))
  #check whether run successfully with the same seed
  set.seed(888)
  mat1 <- aipw$fit()$summary()$result
  set.seed(888)
  mat2 <- aipw$fit()$summary()$result
  expect_output(print(mat1),regexp = "Estimate")
  expect_identical(mat1,mat2)
})
