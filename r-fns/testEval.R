## try doing emu and then taking evals of the cov matrix
source("EmuRbind.R")
source("plotCovReals.R")
library("Matrix")

# npts
npts <- 8
model <- demoModel(npts, lhs=0, rangeMin=0.0, rangeMax=1.5)

setEmulatorOptions(0, 1.9 )
thetas <- callEstimate(model, npts)

results <- callEmulate(model, thetas, npts, rangemin=0.0, rangemax=1.5)

# this is going to be our analytic model for plotting against
# not transparent but the fn is def'd in testRbind
sequence <- seq(0.0, 1.0, length=100)
actual <- data.frame(x=sequence, y=yM(sequence))


# we can plot the results, and it works!
plotResultsTest(model, results, actual, "our model")


## now make the cmatrix
cM <- makeCMatrix(npts, thetas, model$xmodel)
eigs <- eigen(cM)

# this is  just useful for manipulation
vecs <- eigs$vectors
evals <- eigs$values
# note that cM == vecs %*% f %*% t(vecs)
#f <- Matrix(0, npts, npts) <- diag(evals)

## now try something silly
x <- seq(0, 1.5, length.out=200)
for(i in 1:npts){
  y <- rep(NA, 200)
  for(j in 1:200){
    y[j] <- dnorm(x[j], mean=model$xmodel[i], sd=evals[i])
  }
  lines(x,y)
}


