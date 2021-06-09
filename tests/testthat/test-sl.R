#' @title Testing sl.fit and sl.predict wrapper: SuperLearner
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/06/02
test_that("sl.fit & sl.predict: SuperLearner", {
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
  #sl.fit function
  SL_fit <- aipw$sl.fit(X=aipw$.__enclos_env__$private$Q.set,
                        Y=aipw$.__enclos_env__$private$Y,
                        SL.library=sl.lib,
                        CV=list())
  expect_identical(class(SL_fit),"SuperLearner")
  #sl.pred function
  SL_pred <- aipw$sl.predict(SL_fit, newdata = aipw$.__enclos_env__$private$Q.set)
  expect_true(is.numeric(SL_pred))
  expect_equal(length(SL_pred),100)
})

