test_that("aipw input data", {
  expect_error(
    aipw_input(Y=rep(1,100),
               A=rep(1,80),
               W=rep(1,100),
               Q.SL.library=c("SL.mean","SL.glm"),
               g.SL.library=c("SL.mean","SL.glm"),
               k_split = 5,verbose = TRUE)
    )

})

