#' Providing parameter estimates from sensitivity analysis function
#'
#' `ghost_par_ests()` Selects certain parameter estimates from the output of the sensitivity analysis.
#' @param step3 The object returned from SA_step3.
#' @param parameter_label The label used for the parameter in the lavaan code.
#' @param remove_NA Remove rows for combinations of phantom variable parameters that resulted in inadmissable solutions in lavaan.
#' @returns A dataframe of the parameter estimates from the lavaan model.
#' @export


ghost_par_ests <- function(step3,
                           parameter_label,
                           remove_NA = FALSE) {
  resultlist <- step3[[1]]
  combos <- step3[[4]]
  parmatrix <- matrix(NA, nrow = 0, ncol = 10)
  comboindex <- matrix(NA, nrow = 0, ncol = 1)

  for (i in 1:nrow(combos)) {
    reptemp = i
    if (length(resultlist[[i]]) > 1) {
      pars <- resultlist[[i]]
      partemp <- pars[which(pars$"label" == parameter_label),]
    }
    else {
      partemp <- rep(NA, 10)
    }
    parmatrix <- rbind(parmatrix, partemp)
    comboindex <- rbind(comboindex, reptemp)
  }

  out <- cbind(comboindex, parmatrix)
  colnames(out) <-
    c("rep",
      "lhs",
      "op",
      "rhs",
      "label",
      "est",
      "se",
      "z",
      "pvalue",
      "ci.lower",
      "ci.upper")
  out <- cbind(combos, out)

  if (remove_NA == TRUE) {
    out <- na.omit(out)

  }

  return(out)
}


