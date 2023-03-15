#' Sensitivity Analysis Function Step 1
#'
#' `SA_step1()` identifies the phantom variables and generates names for their covariance parameters that will be used in SA_step2().
#' @param lavoutput The lavaan output object output from lavaan functions sem() or lavaan() when fitting your observed model.
#' @param mod_obs A lavaan syntax for the observed model.
#' @param mod_phant A lavaan syntax for the phantom variable model.
#' @returns a list containing the names of all phantom covariance parameters.
#' @export
#'

SA_step1 <- function(lavoutput, #lavaan object
                     mod_obs, # lavaan syntax for observe model
                     mod_phant #model with phantom variables
                     ){

fit = lavoutput

C <- fit@SampleStats@cov[[1]]
colnames(C) <- fit@Data@ov.names[[1]]
var_obs <- diag(C)

obsvar<- unique(c(lavaan::lavParseModelString(mod_obs)$lhs,lavaan::lavParseModelString(mod_obs)$rhs))
SAvar <- unique(c(lavaan::lavParseModelString(mod_phant)$lhs,lavaan::lavParseModelString(mod_phant)$rhs))


# find the variables that are new in the phantom model -- those are the phantom vars
phantom_names <- setdiff(SAvar,obsvar)

newmat <- matrix(NA,nrow=length(SAvar),ncol=length(SAvar))
newmat[1:length(obsvar),1:length(obsvar)]= cov2cor(C)

#newnames <- c(fit@Data@ov.names[[1]],phantom_names)
#simpnames <- gsub(mediator,"M",oldnames[[2]])

colnames(newmat)=c(fit@Data@ov.names[[1]],phantom_names)
rownames(newmat)=colnames(newmat)

var_phant <- c(var_obs,c(rep(1,length(phantom_names))))

matrix = newmat

# put cov names into the matrix
parname <- c()
for (i in 1:nrow(newmat)) {
  for (j in 1:ncol(newmat)){
    if (is.na(newmat[i,j]) & i==j) {
      tn <-(paste0("Var",rownames(newmat)[i],colnames(newmat)[j]))
      newmat[i,j]=tn
      parname <- c(parname,tn)
    }
    else if (is.na(newmat[i,j]) & i>j) {
      tn = (paste0("Cov",rownames(newmat)[i],colnames(newmat)[j]))
      newmat[i,j]=tn
      parname <- c(parname,tn)
    }
    else if (is.na(newmat[i,j]) & i<j) {
     tn = (paste0("Cov",rownames(newmat)[j],colnames(newmat)[i]))
     newmat[i,j]=tn
     parname <- c(parname,tn)
    }
  }
}

namemat <- newmat
cov_names <- colnames(newmat)
namemat <- outer(cov_names, cov_names, function(x, y) {
  ifelse(x == y, paste0("Var", x),
         ifelse(x < y, paste0("Cov", x, y),
                paste0("Cov", x, y)))
})

obsname <- setdiff(unique(parname), namemat[lower.tri(namemat)])

message(paste("Here are the phantom covariance matrix parameters:"))
print(paste0(sort(unique(parname)),collapse='","'),quote=FALSE)
message("Choose the names of the phantom covariances that you want to fix to single values and put in a vector. These will be used for the fixed_names argument in the SA_step2 function.  The phantom covariance parameters that you want to vary should be put in a list and used as the test_names argument. ")
message("Here are the observed covariance matrix parameters:")
#print(sort(namemat[lower.tri(namemat)]))
print(sort(setdiff(namemat[lower.tri(namemat)],unique(parname))))
message("Choose which values you want to use for your fixed parameters and put their names in a vector (fixed_values). Make sure the order is the same for both vectors.")
return(list(matrix,parname,namemat,newmat,mod_phant,var_phant))
}





