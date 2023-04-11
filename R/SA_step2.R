#' Step 2 of sensitivity analysis function
#'
#' `SA_step2()` is used to assign values to the phantom covariances. There are three options for assigning values to the phantom covariances: 1. Fix phantom covariance to a numeric value (i.e., 0 or 1), 2. Fix phantom covariance to be equal to another covariance, or 3. Test different values for the phantom covariance.
#' @param fixed_names this is a vector of the covariance parameter names that will be fixed to a single value
#' @param fixed_values A vector containing either single values or the names of other known parameters that you will fix each parameter in fixed_names to.
#' @param test_names A list of vectors containing the names of the covariances that will be varied. If you wish to constrain certain parameters to be equal, you will need to put the names of these parameters in the same vector.
#' @param test_values A list of seq() to try for each parameter.  If you are constraining certain parameters to be equal, you will put these in a list. See example.
#' @param step1 The object created in `SA_step1()`
#' @returns A list containing test covariance matrices that the phantom model will be fit to.
#' @import corpcor
#' @import lavaan
#' @export
SA_step2 <- function(fixed_names, #covariances fixed to single values
                     fixed_values, #values that you fix variables to
                     test_names=NULL, # list of covariances that will be varied. put parameters in same list if you want them to be equal
                     test_values=NULL, # list of values to try for each parameter
                     step1 # previous step
) {
  fixed=fixed_names
  ref=fixed_values
  matrix=step1[[1]]
  parname = step1[[2]]
  namemat = step1[[3]]
  newmat = step1[[4]]
  mod_phant = step1[[5]]
  var_phant = step1[[6]]

  # reference names indices
  ind_ref <- sapply(fixed_values, function(x) {which(namemat == x,arr.ind=TRUE)})
  vals <- c(rep(0,length(ind_ref)))


  for (i in 1:length(ind_ref)){
    #print(length(ind[[i]]))
    # if(length(ind[[i]])<1)
    if(is.na(ind_ref[[i]][2]))
    {print(is.na(ind_ref[[i]][2]))
      vals[i]<-as.numeric(fixed_values[i])
      print(vals[i])}
    else {vals[i] <-(matrix[ind_ref[[i]]]) }
  }

  # put reference values into covariance matrix
  for (i in 1:length(fixed)){
    index = which(newmat==fixed[i],arr.ind=TRUE)[1,]
    matrix[index[[1]],index[[2]]] = as.numeric(vals[i])
    matrix[index[[2]],index[[1]]] = as.numeric(vals[i])
    # matrix <- sub(fixed[i],as.numeric(vals[i]),newmat)
  }

  # variables that are NA
  naind <- which((is.na(matrix)&lower.tri(matrix)),arr.ind=TRUE)

  #names of remaining parameters that don't have values
  name_na <- namemat[naind]

  # check if there are any remaining NA values

  for (j in 1:nrow(naind)){
    print(name_na[j] %in% unlist(test_names))
  }


  matrix[which((is.na(matrix)&lower.tri(matrix)))]

  # if no variables with custom ranges are entered
  if (is.null(test_values) & is.null(test_names)) {
    saparname <- newmat[naind]
    combocols <- nrow(naind)
    range <- seq(-.3,.3,.1)
    combos <- eval(parse(text= paste("crossing(",paste(rep("range",combocols),collapse=","),")")))
    colnames(combos) <- saparname
    combos <- as.data.frame(combos)

    corlist <-
      rep(list((list(
        matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
      ))), nrow(combos))

    for (i in 1:nrow(combos)){
      tmat <- matrix
      tmat[naind]=unlist(combos[i,])
      tmat[upper.tri(tmat)]<-t(tmat)[upper.tri(tmat)]
      corlist[[i]][[1]] = tmat
      corlist[[i]][[2]] = corpcor::is.positive.definite(tmat)
    }

    covlist <-
      rep(list((list(
        matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
      ))), nrow(combos))

  } else if (!is.null(test_values) & !is.null(test_names)) {
    #combos <- reduce(test_values,crossing)
    combos <- expand.grid(test_values)
    #colnames(combos) <- sapply(test_names,"[[",1) # pick first element of each entry in names list
    if (length(test_names)==1){ #if there is one vector of test names then the columns need to correspond to those params
      colnames(combos) <- unlist(test_names)
    }
    else {
    colnames(combos) <- lapply(test_names,paste0,collapse=",") # combines names of params
    }
    combos <- as.data.frame(combos)
    corlist <-
      rep(list((list(
        matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
      ))), nrow(combos))



    for (i in 1:nrow(combos)){
      tmat=matrix
      if (length(test_names)==1){
      tn <- unlist(test_names)
        for (j in 1:length(tn)){
          tmat[(which(namemat==tn[j],arr.ind=TRUE))]=combos[i,j]
        }
    } else {
      for (j in 1:length(test_names)){
        unlist(test_names)
        for (k in 1:length(test_names[[j]])){
          tmat[(which(namemat==test_names[[j]][k], arr.ind=TRUE))]=combos[i,j]
        }
      }
    }
      tmat[upper.tri(tmat)]<-t(tmat)[upper.tri(tmat)]
      corlist[[i]][[1]] = tmat
      corlist[[i]][[2]] = corpcor::is.positive.definite(tmat)
    }

  } else if (!is.null(test_values) & is.null(test_names)) { message("Error: You must provide lists for BOTH (test_values) and (test_names), or leave these null.  You have provided testvalues but have not specified parameters in test_names. ")
  } else if (is.null(test_values) & !is.null(test_names)) {message("Error: You must provide lists for BOTH (test_values) and (test_names), or leave these null. You have provided parameter names (test_names) but have not specified the values you want to test (test_values).")}




  covlist <-
    rep(list((list(
      matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
    ))), nrow(combos))

  for (i in 1:nrow(combos)){
    covlist[[i]] <- lavaan::cor2cov(R=corlist[[i]][[1]], sd=sqrt(var_phant))
  }

  return(list(mod_phant = mod_phant,
              var_phant=var_phant,
              combos=combos,
              correlation_matrix_list = corlist,
              covariance_matrix_list = covlist))

}
