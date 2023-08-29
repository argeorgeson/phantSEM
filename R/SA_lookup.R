#' Lookup Table for Sensitivity Analysis
#'
#' `SA_lookup()` is used to look up the sensitivity analysis results for a two-wave mediation model when provided with the cross-sectional correlations.
#' @param CorXM The observed correlation between predictor X and mediator M.
#' @param CorXY The observed correlation between predictor X and outcome Y.
#' @param CorMY The observed correlation between mediator M and outcome Y.
#' @import dplyr
#' @returns Results of a sensitivity analysis with varying cross-lagged and autoregressive correlations.
#' @export
#' @examples
#' # specify correlations
#' xm <- .2
#' xy <- .3
#' my <- .4
#'
#' output <- SA_lookup(
#'   CorXM = xm,
#'   CorXY = xy,
#'   CorMY = my
#' )
SA_lookup <- function(
    CorXM,
    CorXY,
    CorMY) {
  # we have to round these -- may change in future if lookup table is larger
  xm <- round(CorXM, 1)
  xy <- round(CorXY, 1)
  my <- round(CorMY, 1)


  obs2 <- obs2
  cond_id <- obs2 %>% filter(round(obs2$corxm, 1) == xm, round(obs2$corxy, 1) == xy, round(obs2$cormy, 1) == my)

  # check to see if
  if (nrow(cond_id) == 0) {
    stop("The correlation matrix you have input results in a non-positive definite covariance matrix. This may be due to rounding. Please use `SA_step1`, `SA_step2` and `SA_step3")
  }


  filt_results <- results_noNA[which(results_noNA$cond == cond_id[, 4]), ]

  return(filt_results)
}
