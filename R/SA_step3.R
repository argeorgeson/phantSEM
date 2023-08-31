#' Step 3 of sensitivity analysis function
#'
#' `SA_step3()` computes the parameter estimates in your phantom model defined in step 1 for the different values provided.
#' @param step2 The object returned from SA_step2.
#' @param n The sample size.
#' @returns A list of parameter estimates from each test covariance matrix.
#' @import lavaan
#' @importFrom lavaan sem
#' @export
#' @examples
#' #' @examples
# covariance matrix
#' covmatrix <- matrix(c(
#'   0.25, 0.95, 0.43,
#'   0.95, 8.87, 2.66,
#'   0.43, 2.66, 10.86
#' ), nrow = 3, byrow = TRUE)
#' colnames(covmatrix) <- c("X", "M2", "Y2")
#'
#' # lavann syntax for observed model
#' observed <- " M2 ~ X
#'              Y2 ~ M2+X "
#'
#' # lavaan output
#' obs_output <- lavaan::sem(model = observed, sample.cov = covmatrix, sample.nobs = 200)
#'
#' # lavaan syntax for phantom variable model
#' phantom <- " M2 ~ M1 + Y1 + a*X
#'                Y2 ~ M1 + Y1 + b*M2 + cp*X "
#'
#' Step1 <- SA_step1(
#'   lavoutput = obs_output,
#'   mod_obs = observed,
#'   mod_phant = phantom
#' )
#'
#' phantom_assignment <- list(
#'   "CovM1X" = 0,
#'   "CovY1M1" = "CovY2M2",
#'   "CovY1X" = 0,
#'   "VarM1" = 1,
#'   "VarY1" = 1,
#'   "CovM1M2" = seq(.4, .6, .1),
#'   "CovY1Y2" = "CovM1M2",
#'   "CovY1M2" = seq(.1, .3, .1),
#'   "CovM1Y2" = "CovY1M2"
#' )
#' Step2 <- SA_step2(
#'   phantom_assignment = phantom_assignment,
#'   step1 = Step1
#' )
#' Step3 <- SA_step3(
#'   step2 = Step2,
#'   n = 200
#' )
#'
SA_step3 <- function(step2, n) {
  combos <- step2[[3]]
  mod_phant <- step2[[1]]
  var_phant <- step2[[2]]
  corlist <- step2[[4]]
  covlist <- step2[[5]]

  # indices of NA in matrix

  # naind <- which((is.na(matrix)&lower.tri(matrix)),arr.ind=TRUE)



  sumlist <- list(NA)
  for (i in 1:nrow(combos)) {
    if (corlist[[i]][2] == TRUE) {
      sumlist[[i]] <- lavaan::parameterEstimates(lavaan::sem(model = mod_phant, sample.cov = covlist[[i]], sample.nobs = n))
    } else {
      sumlist[[i]] <- c("NPD")
    }
  }
  return(list(sumlist, corlist, covlist, combos))
}
