# the power-exp covariance function with
# the power(alpha) set to 1.9
# note that the nugget thetas[3] gets added in when we
# make the covariance matrix
covFn <- function(z1,z2,thetas){
  beta <- thetas[4]^(1.9)
  cov <- exp((-1.0/2.0)*(abs(z1-z2))^(1.9)/beta)
  ans <- thetas[1]*cov#+thetas[2]  
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
