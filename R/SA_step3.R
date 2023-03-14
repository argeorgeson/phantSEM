#' Step 3 of sensitivity analysis function 
#' 
#' This function computes the parameter estimates in your phantom model defined in step 1 for the different values provided. 
#' @param step2 the object returned from SA_step2. 
#' @param n the sample size. 
#' SA_step3()


SA_step3 <- function(step2,n){
  combos = step2[[3]]
  mod_phant=step2[[1]]
  var_phant=step2[[2]]
  corlist = step2[[4]]
  covlist = step2[[5]]
  
  # indices of NA in matrix 
  
  naind <- which((is.na(matrix)&lower.tri(matrix)),arr.ind=TRUE)
  
  
  
  sumlist <- list(NA)
  for (i in 1:nrow(combos)){
    if (corlist[[i]][2]==TRUE){
      sumlist[[i]]=parameterEstimates(sem(model=mod_phant,sample.cov = covlist[[i]],sample.nobs = n))
    } else {sumlist[[i]] = c("NPD")}
  }
  return(list(sumlist,corlist,covlist, combos))
}