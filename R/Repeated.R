# Create a decorator that allows AIPW class repeated cross-fitting
Repeated <- R6Class("RepeatedFit",
                          inherit = AIPW,
                          public = list(
                            aipw_obj = NULL,
                            num_reps = NULL,
                            stratified_fitted = NULL,
                            repeated_estimates = NULL,
                            repeated_results = NULL,
                            result = NULL,

                            initialize = function(aipw_obj) {
                              self$aipw_obj = aipw_obj
                              private$verbose = self$aipw_obj$.__enclos_env__$private$verbose
                              private$Y.type = self$aipw_obj$.__enclos_env__$private$Y.type
                              self$aipw_obj$.__enclos_env__$private$verbose = FALSE

                              #-------check if future.apply is loaded otherwise lapply would be used.------#
                              if (any(names(sessionInfo()$otherPkgs) %in% c("future.apply"))){
                                private$use.f_lapply = TRUE
                              } else {
                                private$use.f_lapply = FALSE
                              }
                            },

                            repfit = function(num_reps, stratified){
                              self$num_reps = num_reps
                              self$stratified_fitted = stratified

                              iter = 1:self$num_reps
                              self$repeated_estimates =
                                private$.f_lapply(iter,function(i,...){
                                  self$aipw_obj$fit()$summary()
                                  estimates_count = 3
                                  Estimand_label = c("Risk of exposure", "Risk of control","Risk Difference")

                                  if (private$Y.type == 'binomial'){
                                    estimates_count = estimates_count +2
                                    Estimand_label = c(Estimand_label,"Risk Ratio", "Odds Ratio")
                                  }

                                  if (self$stratified_fitted) {
                                    estimates_count = estimates_count +2
                                    Estimand_label = c(Estimand_label, "ATT Risk Difference","ATC Risk Difference")
                                  }

                                  estimates = data.frame(do.call(rbind,self$aipw_obj$estimates[1:estimates_count]))
                                  estimates$Estimand = rownames(estimates)
                                  estimates$Estimand = factor(estimates$Estimand,
                                                              levels = estimates$Estimand,
                                                              labels = Estimand_label)
                                  return(estimates)
                                })
                              names(self$repeated_estimates) = iter
                              self$repeated_estimates = data.frame(do.call(rbind, self$repeated_estimates))
                              colnames(self$repeated_estimates) = c("Estimate","SE","Estimand")

                            },

                            summary = function(){
                              self$repeated_results = split(self$repeated_estimates, self$repeated_estimates$Estimand)
                              self$repeated_results = lapply(self$repeated_results,
                                                             function(x) {
                                                               if(unique(x$Estimand) %in% c("RR","OR")){
                                                                 x$Estimate = log(x$Estimate)
                                                                 median_adjusted = private$get_median_variance(x$Estimate,x$SE)
                                                                 median_adjusted[
                                                                   names(median_adjusted) %in%
                                                                     c("Median Estimate","95% LCL Median Estimate", "95% UCL Median Estimate")
                                                                   ] = exp(median_adjusted[
                                                                     names(median_adjusted) %in%
                                                                       c("Median Estimate","95% LCL Median Estimate", "95% UCL Median Estimate")
                                                                   ] )
                                                               } else{
                                                                 median_adjusted = private$get_median_variance(x$Estimate,x$SE)
                                                               }

                                                               return(median_adjusted)
                                                             })
                              self$result = do.call(rbind, self$repeated_results)

                              if (private$verbose){
                                print(self$result,digit=3)
                              }

                            }
                            ),

                       private = list(
                         get_median_variance = function(est, se){
                           est_median = median(est)
                           est_var = se^2
                           se_median = median(se)
                           est_median_var = median( est_var +  (est - est_median)^2)
                           est_median_se = sqrt(est_median_var)
                           ci_median_se = get_ci(est_median,est_median_se)
                           res = c(est_median,se_median,est_median_se,ci_median_se)
                           names(res) = c("Median Estimate", "Median SE", "SE of Median Estimate","95% LCL Median Estimate","95% UCL Median Estimate")

                           return(res)
                         }
                       )
)

