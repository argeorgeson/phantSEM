#' Step 2 of sensitivity analysis function
#'
#' `SA_step2()` is used to assign values to the phantom covariances. There are three options for assigning values to the phantom covariances: 1. Fix phantom covariance to a numeric value (i.e., 0 or 1), 2. Fix phantom covariance to be equal to another covariance, or 3. Test different values for the phantom covariance.
#' @param phantom_assignment A list of all phantom parameter names (copied from `SA_step1()` output) which assigns them to be equal to ONE of the following: 1) an observed parameter name, 2) a single numeric value, 3) a sequence of values, or 4) another phantom variable that has been set equal to 1-3.
#' @param step1 The object created in `SA_step1()`
#' @returns A list containing test covariance matrices that the phantom model will be fit to.
#' @import corpcor
#' @import lavaan
#' @export
SA_step2 <- function(phantom_assignment, # list of all phantom parameter names which tells the function what they should be equal to
                     step1 # previous step
) {

  pa = phantom_assignment
  matrix_template=step1[[1]]
  parname = step1[[2]]
  namemat = step1[[3]]
  newmat = step1[[4]]
  mod_phant = step1[[5]]
  var_phant = step1[[6]]
  cov_map = step1[[7]]

  # covariance names fixed to single numeric value
  fvaltable <- data.frame(matrix(ncol = 2,nrow=0))
  for (i in seq_along(pa)) {
    if (is.numeric(pa[[i]]) && length(pa[[i]]) == 1) {
      fvaltable <- rbind(fvaltable, c(names(pa)[i], pa[[i]]))
    }
  }
  colnames(fvaltable) <- c("name", "value")



# covariance names set to another named value
  charvaltable <- unlist(pa[sapply(pa,is.character)])
  charvaltable <-setNames(charvaltable,NULL)


  #this returns the rows that have an observed parameter value
char_obs <- subset(cov_map,covname %in% charvaltable & !is.na(val))
 # covariances equal to other values
eqseq <- subset(cov_map,covname %in% charvaltable & is.na(val))[,1]

# covariances that are going to be varied
seqvaltable <- lapply(pa, function(x) if (is.numeric(x) & length(x)>1) x )
seqvaltable <- seqvaltable[!sapply(seqvaltable,is.null)]


# covariances that are fixed to observed covariances -- by name
fnametable <- data.frame(matrix(ncol = 2,nrow=0))
for (i in seq_along(pa)) {
  if (is.character(pa[[i]]) && length(pa[[i]]) == 1 && !(pa[[i]] %in% eqseq)) {
    fnametable <- rbind(fnametable, c(names(pa)[i], pa[[i]]))
  }
}
colnames(fnametable) <- c("name", "value")

# this is the fixed name and values
fixed <- rbind(fnametable,fvaltable)

# covariances that are fixed to other phantom covariances
phantomnametable <- data.frame(matrix(ncol = 2,nrow=0))
for (i in seq_along(pa)) {
  if (is.character(pa[[i]]) && length(pa[[i]]) == 1 && (pa[[i]] %in% eqseq)) {
    phantomnametable <- rbind(phantomnametable, c(names(pa)[i], pa[[i]]))
  }
}
colnames(phantomnametable) <- c("name","value")




# this creates list of parameter names and values
testnamelist <- list(NA)
testnametable <- data.frame(matrix(nrow=0,ncol=1))
#counter
j=0
for (i in seq_along(pa)) {
  if (is.numeric(pa[[i]]) && length(pa[[i]]) > 1 ) {
    j=j+1
    print(j)
    testnametable <- rbind(testnametable, names(pa)[i])
    testnamelist[[j]] <- list(name=names(pa)[i],values=pa[[i]])
  }
}
colnames(testnametable) <- c("value")


# make sure the phantom variable name is correct
for (i in 1:nrow(phantomnametable)) {
  phantomnametable[i,2] %in% testnametable
}

#make sure that when there are phantom variables equal to other phantom variables, the phantom variable name is correct
stopifnot("The phantom variable name that another phantom variable was set equal to is not spelled correctly. Please compare to output from SA_step1" =
         all(phantomnametable[,2] %in% testnametable))


# if all phantom variables are fixed to numeric values
if(length(unique(parname))==nrow(fvaltable)) {


} else {

# Create a new list
new_testnamelist <- list()
# Loop over the list of lists
for (i in seq_along(testnamelist)) {

  # Get the name and values of the current list
  name <- testnamelist[[i]]$name
  values <- testnamelist[[i]]$values

  # Find the rows in the data frame that matches the value
  rows <- phantomnametable[phantomnametable$value == name, ]

  # Add the name(s) from the data frame to the current list
  new_testnamelist[[i]] <- list(name = c(name,rows$name), values = values)
}

}
#### assigning names that were previously function arguments ####
fixed_names=fixed[,1]
fixed_values=fixed[,2]
test_names= lapply(new_testnamelist, '[[',1)
test_values <- lapply(new_testnamelist, '[[',2)
fixed=fixed_names
ref=fixed_values

  # reference name indices
  ind_ref <- sapply(fixed_values, function(x) {which(namemat == x,arr.ind=TRUE)})
  vals <- c(rep(0,length(ind_ref)))


  for (i in 1:length(ind_ref)){
    #print(length(ind[[i]]))
    # if(length(ind[[i]])<1)
    if(is.na(ind_ref[[i]][2]))
    {print(is.na(ind_ref[[i]][2]))
      vals[i]<-as.numeric(fixed_values[i])
      print(vals[i])}
    else {vals[i] <-(matrix_template[ind_ref[[i]]]) }
  }

  # put reference values into covariance matrix
  for (i in 1:length(fixed)){
    index = which(newmat==fixed[i],arr.ind=TRUE)[1,]
    matrix_template[index[[1]],index[[2]]] = as.numeric(vals[i])
    matrix_template[index[[2]],index[[1]]] = as.numeric(vals[i])
    # matrix <- sub(fixed[i],as.numeric(vals[i]),newmat)
  }

  # variables that are NA
  naind <- which((is.na(matrix_template)&lower.tri(matrix_template)),arr.ind=TRUE)

  #names of remaining parameters that don't have values
  name_na <- namemat[naind]

  # check if there are any remaining NA values
  for (j in 1:nrow(naind)){
    print(name_na[j] %in% unlist(test_names))
  }

  # if fixed values are given for all of the phantom variables
  if(length(unique(parname))==nrow(fvaltable)) {

    corlist<-
      rep(list((list(
        matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
      ))), 1)

    corlist[[1]] <- matrix_template

    covlist <-
      rep(list((list(
        matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
      ))), 1)

      covlist[[1]] <- lavaan::cor2cov(R=corlist[[1]], sd=sqrt(var_phant))

  } else {

  # if no variables with custom ranges are included on the list -- default is to use range -.3 to .3
 if (length(naind)>0){
  # if (is.null(test_values) & is.null(test_names)) {
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
      tmat <- matrix_template
      tmat[naind]=unlist(combos[i,])
      tmat[upper.tri(tmat)]<-t(tmat)[upper.tri(tmat)]
      corlist[[i]][[1]] = tmat
      corlist[[i]][[2]] = corpcor::is.positive.definite(tmat)
    }

    covlist <-
      rep(list((list(
        matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
      ))), nrow(combos))

    # if there are test values given
  } else if (length(test_values)>0 & length(test_names)>0) {
    #combos <- reduce(test_values,crossing)
    combos <- expand.grid(test_values)

    ## this code is causing issues
    #if (length(test_names)==1){ #if there is one vector of test names then the columns need to correspond to those params
    #  colnames(combos) <- unlist(test_names)
    #}
    ##

    #else { #this code is causing issues
    colnames(combos) <- lapply(test_names,paste0,collapse=",") # combines names of params
    #} #this code is causing issues
    combos <- as.data.frame(combos)
    corlist <-
      rep(list((list(
        matrix(NA, nrow = nrow(newmat), ncol = ncol(newmat)), c(NA)
      ))), nrow(combos))



    for (i in 1:nrow(combos)){
      tmat=matrix_template
      if (length(test_names)==1){
      tn <- unlist(test_names)
        for (j in 1:length(tn)){
          tmat[(which(namemat==tn[j],arr.ind=TRUE))]=combos[i,] #combos would only be a vector here so do not index column
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
  }

  return(list(mod_phant = mod_phant,
              var_phant=var_phant,
              combos=combos,
              correlation_matrix_list = corlist,
              covariance_matrix_list = covlist))

}
