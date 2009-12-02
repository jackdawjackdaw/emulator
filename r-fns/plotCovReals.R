# the power-exp covariance function with
# the power(alpha) set to 1.9
# note that the nugget thetas[3] gets added in when we
# make the covariance matrix
covFn <- function(z1,z2,thetas){
  beta <- thetas[4]^(1.9)
  cov <- exp((-1.0/2.0)*(abs(z1-z2))^(1.9)/beta)
  ans <- thetas[1]*cov
}

# construct the covariance matrix of a given input model
makeCMatrix <- function(m, thetas, xmodel){
  covM <- diag(m) # (e^0)
  for(j in 1:(m-1)) 
    for(i in (1+j):m) 
      covM[i,j] <- covM[j,i] <- covFn(xmodel[i], xmodel[j], thetas)
  # add in the nugget
  diag(covM) <- 1 + thetas[3]
  covM
}

## takes into account the conditioning on the covariance
## due to the model points (not just the hyper params)
## npts -> number of points in bigmodel
## m -> number of pionts in xmodel
## thetas -> hyperparams of the covfn (shared by bigmodel and xmodel)
## xmodel -> the model data used by the emulator
## bigmodel -> the grid of points we want to sample the emulator distribution at
##
## adjustedC = cBig - Kij cModel^-1 t(Kij)
makeAdjustedC <- function(npts, m, thetas, xmodel, bigXmodel){
  ## make the covariance matrix of the original model, and cholesky it straight off
  #browser()
  cModel <- makeCMatrix(m, thetas, xmodel)
  cModelChol <- chol(cModel)
  cModelInv <- chol2inv(cModelChol)
  ## this is the cov matrix for the grid of points we need to
  ## make a nice plot of the samples
  cBig <- makeCMatrix(npts, thetas, bigXmodel)

  ##
  Kij <- matrix(0, npts, m)
  for(i in 1:npts){
    for(j in 1:m){
      ## ahh unfortunate notation here
      Kij[i,j] <- covFn(bigXmodel[i],xmodel[j], thetas)
    }
  }

  cFinal <- matrix(0, npts, npts)
  #browser()
  cFinal <- cBig - (Kij %*% (cModelInv %*% t(Kij)))
  cFinal
}
