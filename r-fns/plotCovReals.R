# the power-exp covariance function with
# the power(alpha) set to 1.9
# note that the nugget thetas[3] gets added in when we
# make the covariance matrix
covFn <- function(z1,z2,thetas, alpha=1.9){
  beta <- thetas[4]^(alpha)
  cov <- exp((-1.0/2.0)*(abs(z1-z2))^(alpha)/beta)
  ans <- thetas[1]*cov
}

# construct the covariance matrix of a given input model
makeCMatrix <- function(m, thetas, xmodel, alpha=1.9){
  covM <- diag(m) # (e^0)
  for(j in 1:(m-1)) 
    for(i in (1+j):m) 
      covM[i,j] <- covM[j,i] <- covFn(xmodel[i], xmodel[j], thetas, alpha)
  # add in the nugget
  diag(covM) <- 1 + thetas[3]
  covM
}


## make some reals of a gp with diff alpha using the
makeReal <- function(alpha, npts, thetaLength=0.1){
  grid <- seq(0,1, length.out=npts)
  thetas <- c(1.0, 0.0, 0.0, thetaLength)
  covMatrix <- makeCMatrix(npts, thetas, grid, alpha)
  A <- chol(covMatrix)
  z <- rnorm(npts)
  ans <- data.frame(x=grid,y=t(A)%*%z)
}


## this one is less informative
plotRealsLength <- function(){
  a <- makeReal(1.9, 100, 0.1)
  b <- makeReal(1.9, 100, 0.05)
  c <- makeReal(1.9, 100, 0.3)
  plot(a, col="black", ylim=range(-2.0,2.0), type="l", lwd=2)
  lines(b, col="blue", lwd=2)
  lines(c, col="green", lwd=2)
  title(main="Realisations of a Gaussian Process")
  legend(x=0.7, y=2, bg="white", legend=c('theta_2=0.1', 'theta_2=0.05', 'theta_2=0.3'),
       col=c('black', 'blue', 'green'),
       lwd=2)
  grid()
}  

plotReals <- function(){
  a <- makeReal(1.99, 100)
  b <- makeReal(1.7, 100)
  c <- makeReal(1.0, 100)
  plot(a, col="black", ylim=range(-2.0,2.0), type="l", lwd=2)
  lines(b, col="blue", lwd=2)
  lines(c, col="green", lwd=2)
  title(main="Realisations of a Gaussian Process")
  legend(x=0.7, y=2, bg="white", legend=c('alpha=1.99', 'alpha=1.7', 'alpha=1.0'),
       col=c('black', 'blue', 'green'),
       lwd=2)
  grid()

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
