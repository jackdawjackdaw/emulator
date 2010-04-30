library("Matrix")
# try and implement the sa functions
# these are derived with the idea that regression kernel is a squared-exponential
# with a diagonal covariance structure Sigma
# and a simple linear regression vector h(x) = (1, x)

# the one point R vector
# p is the index of the parameter we're doing the SA over
# x is the vector of param values
# nq is the number of regression functions, this sets the length of Rp
rpOne <- function(x, p, nq){
  ans <- rep(0, nq)
  ans[1] <- 1
  ans[p+1] <- x
  ans

}

# the trivial case of no set
rZero <- function(nq){
  ans <- rep(0, nq)
  ans[1] <- 1
  ans
}

# the one point T vector
# x,p as above
# desginMatrix is the matrix of of design points
# ndes is number of dimensions in the design, length of a design vector
# paramVariances is the diagonal matrix of variances of our inputs (set apriori)
# gpEmuVariances, the diagonal matrix of variances from the emulator (thetas)
tpOne <- function(x, p, designMatrix, nmodelpts, ndim, paramVariances, gpEmuVariances){
  exponentFirst <- rep(0, nmodelpts)
  exponentSec <- rep(0, nmodelpts)
  designVecSubP <- rep(0, (ndim-1))
  jVec <- rep(0, (ndim-1))
  gpEmuVariancesSubP <- matrix(0, nrow=ndim-1, ncol=ndim-1) # this is bhat
  paramVariancesSubP <- matrix(0, nrow=ndim-1, ncol=ndim-1) # and this is M

  gpEmuVariancesSubP <- gpEmuVariances[-p, -p]
  paramVariancesSubP <- paramVariances[-p, -p]
  
  for(i in 1:nmodelpts){
    # compared to the notes we have
    # designVec = \tilde{x_i}
    # gpEmuvariances = B
    # gpEmuvariances[p,p] = B_{pp}
    # x[p] = xp
    exponentFirst[i] = -(1.0/2.0)*(designMatrix[i,] %*% gpEmuVariances %*% designMatrix[i,]  + x*gpEmuVariances[p,p]*x - 2*designMatrix[i,p]*gpEmuVariances[p,p]*x)
  }

  detBM <- det(gpEmuVariancesSubP + paramVariancesSubP)
  # calling solve without another arg will default to
  # inverting the matrix a  
  invBM <- solve(gpEmuVariancesSubP + paramVariancesSubP)
  sqrtVal <- sqrt(det(paramVariancesSubP)/(detBM))

  # now do the J vector part
  # jvec = 2* \tilde{x_i_{-p}} * \hat{B}
  for(i in 1:nmodelpts){
    # create \tilde{x_i_{-p}}
    # need to do this for each 
    designVecSubP <- designMatrix[i, -p]
    # jVec = \tilde{x_i_{-p}} * Beta
    jVec <- designVecSubP %*% gpEmuVariancesSubP
    # not sure about the transpose
    exponentSec[i] <- (-1.0/2.0)*(jVec %*% invBM %*% t(jVec))
  }
  # put it all together
  # this will be a vector of length nmodelpts
  answer <- sqrtVal*exp(exponentFirst)*exp(exponentSec)
}

tZero <- function(designMatrix, nmodelpts, ndim, paramVariances, gpEmuVariances){
  jVec <- rep(0, ndim)
  exponentFirst <- rep(0, nmodelpts)
  exponentSec <- rep(0, nmodelpts)
  for(i in 1:nmodelpts){
    exponentFirst[i] <- -(1.0/2.0)*(designMatrix[i,] %*% gpEmuVariances %*% designMatrix[i,])
  }
  sqrtVal<- sqrt(det(gpEmuVariances) / (det(paramVariances + gpEmuVariances)))
  invBM <- solve(paramVariances + gpEmuVariances)
  for(i in 1:nmodelpts){
    # not sure on the sign here
    jVec <- (designMatrix[i,] %*% gpEmuVariances)
    exponentSec[i] <- -(1.0/2.0)*(jVec %*% invBM %*% t(jVec))
  }

  answer <- sqrtVal*exp(exponentFirst)*exp(exponentSec)
}
  

## and now some support functions
## need these to actually evaluate the conditional expectation value etc
makeHMatrix <- function(designMatrix, ndes, nreg){
  hMatrix <- matrix(0, ncol=nreg, nrow=ndes)
  for(i in 1:ndes)
      hMatrix[i,] <- makeHVector(designMatrix[i,], nreg)

  hMatrix
}

makeHVector <- function(xvalues, nreg){
  result <- rep(0, nreg)
  result[1] <- 1
  for(i in 1:(nreg-1)){
    result[i+1] <- xvalues[i]
  }
  result
}

makeCovMatrix <- function(designMatrix, thetas, ndes){
  covM <- diag(ndes)
  for(j in 1:(ndes-1)) 
    for(i in (1+j):ndes) 
      covM[i,j] <- covM[j,i] <- covFn(designMatrix[i,], designMatrix[j,], thetas)
  covM
}
  

covFn <- function(z1,z2,thetas){
  amp <- thetas[1]
  nugg <- thetas[2]
  beta <- thetas[3:length(thetas)]^2
  exponent <- rep(0, length(beta))
  for(i in 1:length(beta))
    exponent[i] <- (z1[i] - z2[i])^2 / beta[i]
  ## and here's the covariance
  amp*exp(-(0.5)*sum(exponent)) + nugg
}

makeBetaVector <- function(hmatrix, invCovMatrix, training){
  denom <- solve(t(hmatrix) %*% invCovMatrix %*% hmatrix)
  numerator <- t(hmatrix) %*% invCovMatrix %*% training
  denom %*% numerator
}

makeEVector <- function(beta, hmatrix, training, invCovMatrix){
  invCovMatrix  %*% (training - hmatrix %*% beta)
}
