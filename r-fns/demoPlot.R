## a demo
source("EmuRbind.R") ## provides the basic interface to the c fns via libRbind
source("testRbind.R") ## provides some model test functions for plotting etc
library("Matrix")

covFn <- function(z1,z2,thetas){
  beta <- thetas[4]^(1.9)
  cov <- exp((-1.0/2.0)*(abs(z1-z2))^(1.9)/beta)
  ## if(abs(z1-z2) < 0.001){
  ##   print("added nugget")
  ##   nugget <- thetas[3]
  ## } else {
  ##   nugget <- 0
  ## }
  ## theats[2] is what makes the matrix
  ## non pos-def, remove from the c-code!
  ans <- thetas[1]*cov#+thetas[2]  
}

makeCMatrix <- function(m, thetas, xmodel){
  covM <- diag(m) # (e^0)
  for(j in 1:(m-1)) 
    for(i in (1+j):m) 
      covM[i,j] <- covM[j,i] <- covFn(xmodel[i], xmodel[j], thetas)
  # add in the nugget
  diag(covM) <- 1 + thetas[3]
  covM
}


# npts
npts <- 8


# this function in testRbind.R creates some model data
# with the yM wobbly sin(exp) fn over the range [0,1]
# the first arg sets the number of points to distribute in that range
# the second arg sets the if we should use lhs or not (we do here)
ourModel <- demoModel(npts, lhs=0)

# now we set the default values for the covariance function
# we'll need to know this later when we want to rebuild the
# C matrix
#setDefaultOps()
setEmulatorOptions(1,2.0)

# now we'll estimate the thetas for our model
# see EmuRind.R for some info on this one
thetas <- callEstimate(ourModel, npts)

print(thetas)

# now we'll make some emulated data from the thetas
# by default this assumes 1 param, 4 thetas, 100 emupts and
# the range [0, 1]
results <- callEmulate(ourModel, thetas, npts)

# this is going to be our analytic model for plotting against
# not transparent but the fn is def'd in testRbind
sequence <- seq(0.0, 1.0, length=100)
actual <- data.frame(x=sequence, y=yM(sequence))


# we can plot the results, and it works!
plotResultsTest(ourModel, results, actual, "our model")


##
## this is ok but really we should do this at quite a few points
## 



# now make the covMatrix
cM <- makeCMatrix(npts, thetas, ourModel$xmodel)
bigpts<- 100
cMHuge <- makeCMatrix(bigpts, thetas, results$emulatedx)

# make the cholesky decomp
# this doesn't seem to work because the cM matrix is not positive def,
# this is worrying, it really should be
f <- chol(cM)
f2 <- chol(cMHuge)
for(i in 1:10){
  z <- rnorm(npts)
  z1 <- rnorm(bigpts)
# now we have some samples with the right correlation
  samples <- t(f) %*% z
  samples2 <- t(f2) %*% z1
  points(results$emulatedx, results$emulatedy+samples2, col="antiquewhite")

}
