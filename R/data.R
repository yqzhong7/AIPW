#' Simulated Observational Study
#'
#' Datasets were simulated using baseline covariates (sampling with replacement) from the Effects of Aspirin in Gestation and Reproduction (EAGeR) study.
#' Data generating mechanisms were described in our manuscript (Zhong et al. (inpreparation), Am. J. Epidemiol.).
#' True marginal causal effects on risk difference, log risk ratio and log odds ratio scales were attached to the dataset attributes (true_rd, true_logrr,true_logor).
#'
#' @docType data
#'
#' @usage data(eager_sim_obs)
#'
#' @format An object of class data.frame with 200 rows and 8 columns:
#' \describe{
#'   \item{sim_Y}{binary, simulated  outcome which is condition on all other covariates in the dataset}
#'   \item{sim_A}{binary, simulated exposure which is conditon on all other covarites expect sim_Y.}
#'   \item{eligibility}{binary, indicator of the eligibility stratum}
#'   \item{loss_num}{count, number of prior pregnancy losses}
#'   \item{age}{continuous, age in years}
#'   \item{time_try_pregnant}{count, months of conception attempts prior to randomization}
#'   \item{BMI}{continuous, body mass index}
#'   \item{meanAP}{continuous, mean arterial blood pressure}
#' }
#' @references Schisterman, E.F., Silver, R.M., Lesher, L.L., Faraggi, D., Wactawski-Wende, J., Townsend, J.M., Lynch, A.M., Perkins, N.J., Mumford, S.L. and Galai, N., 2014. Preconception low-dose aspirin and pregnancy outcomes: results from the EAGeR randomised trial. The Lancet, 384(9937), pp.29-36.
#' @references Zhong, Y., Naimi, A.I., Kennedy, E.H., (In preparation). AIPW: An R package for Augmented Inverse Probability Weighted Estimation of Average Causal Effects. American Journal of Epidemiology
#' @seealso [eager_sim_rct]
"eager_sim_obs"

#' Simulated Randomized Trial
#'
#' Datasets were simulated using baseline covariates (sampling with replacement) from the Effects of Aspirin in Gestation and Reproduction (EAGeR) study.
#'
#' @docType data
#'
#' @usage data(eager_sim_rct)
#'
#' @format An object of class data.frame with 1228 rows and 8 columns:
#' \describe{
#'   \item{sim_Y}{binary, simulated  outcome which is condition on all other covariates in the dataset}
#'   \item{sim_T}{binary, simulated treatment which is condition on eligibility only.}
#'   \item{eligibility}{binary, indicator of the eligibility stratum}
#'   \item{loss_num}{count, number of prior pregnancy losses}
#'   \item{age}{continuous, age in years}
#'   \item{time_try_pregnant}{count, months of conception attempts prior to randomization}
#'   \item{BMI}{continuous, body mass index}
#'   \item{meanAP}{continuous, mean arterial blood pressure}
#' }
#'
#' @references Schisterman, E.F., Silver, R.M., Lesher, L.L., Faraggi, D., Wactawski-Wende, J., Townsend, J.M., Lynch, A.M., Perkins, N.J., Mumford, S.L. and Galai, N., 2014. Preconception low-dose aspirin and pregnancy outcomes: results from the EAGeR randomised trial. The Lancet, 384(9937), pp.29-36.
#' @references Zhong, Y., Naimi, A.I., Kennedy, E.H., (In preparation). AIPW: An R package for Augmented Inverse Probability Weighted Estimation of Average Causal Effects. American Journal of Epidemiology
#' @seealso [eager_sim_obs]

"eager_sim_rct"
