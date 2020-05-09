#' @title Tesing helper functions
#' @section Last Updated By:
#' Yongqi Zhong
#' @section Last Update Date:
#' 2020/05/09
test_that("AIPW helper functions", {
  #get CI
  expect_equal(AIPW:::get_ci(est=0,se=1,ratio = F),c(-1.96,1.96))
  expect_equal(AIPW:::get_ci(est=1,se=1,ratio = T),exp(c(-1.96,1.96)))
  expect_warning(AIPW:::get_ci(est=-1,se=1,ratio = T))

  #private methods: get_RD, get_RR, get_OR, get_sigma_covar
  vec <- rep(1,100)
  Z_norm = sqrt(length(vec))
  sl.lib <- c("SL.mean","SL.glm")
  aipw <-  AIPW$new(Y=vec,
                    A=vec,
                    W.Q =vec,
                    W.g =vec,
                    Q.SL.library=sl.lib,
                    g.SL.library=sl.lib,
                    k_split = 1,verbose = FALSE)
  #get_RD
  expect_identical(as.numeric(aipw$.__enclos_env__$private$get_RD(aipw_eif0 = vec,
                                                                  aipw_eif1 = vec,
                                                                  Z_norm=Z_norm)),
               rep(0,4))
  #get_sigma_covar
  expect_identical(aipw$.__enclos_env__$private$get_sigma_covar(aipw_eif0 = vec,
                                                                aipw_eif1 = vec),
                   matrix(rep(0,4),ncol=2))
  ####To be done by Ashley
  #get_RR
  # sigma_covar <- matrix(c(0,1,1,0),ncol=2)
  # expect_identical(as.numeric(aipw$.__enclos_env__$private$get_RR(aipw_eif0 = vec,
  #                                                                 aipw_eif1 = vec,
  #                                                                 sigma_covar = sigma_covar,
  #                                                                 Z_norm=Z_norm)),
  #                  rep(0,4))
  #get_OR
  # expect_identical(as.numeric(aipw$.__enclos_env__$private$get_OR(aipw_eif0 = vec,
  #                                                                 aipw_eif1 = vec,
  #                                                                 sigma_covar = sigma_covar,
  #                                                                 Z_norm=Z_norm)),
  #                  rep(0,4))


  }
  )
