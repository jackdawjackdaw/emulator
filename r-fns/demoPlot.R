## a demo
source("EmuRbind.R") ## provides the basic interface to the c fns via libRbind
source("testRbind.R") ## provides some model test functions for plotting etc
source("plotCovReals.R") ## makes a covariance matrix
library("Matrix")



# number of model points
npts <- 8

# number of repeat draws to make from the cov fn
nreps <- 5

# what region to do the emulation in
xmin <- 0.0
xmax <- 1.5

# this function in testRbind.R creates some model data
# with the yM wobbly sin(exp) fn over the range [0,1]
# the first arg sets the number of points to distribute in that range
# the second arg sets the if we should use lhs or not (we do here)
ourModel <- demoModel(npts, lhs=1, rangeMin=xmin, rangeMax=xmax)

# now we set the default values for the covariance function
# we'll need to know this later when we want to rebuild the
# C matrix
setDefaultOps()
#setEmulatorOptions(1,2.0)

# now we'll estimate the thetas for our model
# see EmuRind.R for some info on this one
thetas <- callEstimate(ourModel, npts)

print(thetas)

# now we'll make some emulated data from the thetas
# by default this assumes 1 param, 4 thetas, 100 emupts and
# the range [0, 1]
results <- callEmulate(ourModel, thetas, npts, rangemin=xmin, rangemax=xmax)
bigpts<- 500
bigRes <- callEmulate(ourModel, thetas, npts, nemupts=bigpts, rangemin=xmin, rangemax=xmax)

# this is going to be our analytic model for plotting against
# not transparent but the fn is def'd in testRbind
sequence <- seq(xmin, xmax, length=100)
actual <- data.frame(x=sequence, y=yM(sequence))


# we can plot the results, and it works!
# plot a jpg
jpeg("model-showing-samples.jpg", quality=100, bg="white", res=300, width=8, height=5, units="in")
#plotResultsTest(ourModel, results, actual, "our model")


##
## this is ok but really we should do this at quite a few points
## 

# now make the covMatrix
cM <- makeCMatrix(npts, thetas, ourModel$xmodel)

#cMHuge <- makeCMatrix(bigpts, thetas, results$emulatedx)
cMHugeOrig <- makeCMatrix(bigpts, thetas, bigRes$emulatedx)

cMHuge <- makeAdjustedC(bigpts, npts, thetas, ourModel$xmodel, bigRes$emulatedx)

# make the cholesky decomp,
# now have f2' . f2 = cM) up to numeric errors
f2 <- chol(cMHuge)
f3 <- chol(cMHugeOrig)

# setup the dual plot
par(mfrow=c(1,2))

plot(actual$x, actual$y, type="n", xlim=range(xmin,xmax), ylim=range(0.0,6.0), xlab='x', ylab='y')
grid()
for(i in 1:nreps){
  z1 <- rnorm(bigpts)
# now we have some samples with the right correlation
  samples2 <- t(f2) %*% z1
  # this is a color object, the alpha option makes the points quite transparent
  colTest <- rgb(0, 0, 190, alpha=50, maxColorValue=255)
  points(bigRes$emulatedx, bigRes$emulatedy+samples2, col=colTest, pch=16, type="p", cex=0.75 )

}

                                        # this works in 1d now
points(ourModel$xmodel, ourModel$training, ylim=range(0.0,6.0), xlab="x", ylab="y", pch=19)
title(main='Adjusted Sample Cov Matrix')
lines(actual$x, actual$y, col="black", lwd=2, lty=2)
lines(results$emulatedx, results$emulatedy, col="red", lwd=2)
confidence <- rep(NA, length(results$emulatedvar))
for(i in 1:length(results$emulatedvar))
  confidence[i] <- 0.5*sqrt(results$emulatedvar[i])*1.65585

lines(results$emulatedx, results$emulatedy + confidence, col="red", lty=2)
lines(results$emulatedx, results$emulatedy - confidence, col="red", lty=2)

## add a legend
## legend(x=0.8, y=5, legend=c('model', 'emulator', 'confidence', 'samples calculated Correctly'),
##        col=c('black', 'red', 'red', colTest),
##        lwd=2,
##        lty=c(2,1,2,1))


# and now we'll plot the old one on another graph to comapre?
plot(actual$x, actual$y, type="n", xlim=range(xmin,xmax), ylim=range(0.0,6.0), xlab='x', ylab='y')
grid()
for(i in 1:nreps){
  z1 <- rnorm(bigpts)
# now we have some samples with the right correlation
  samples2 <- t(f3) %*% z1
  # this is a color object, the alpha option makes the points quite transparent
  colTest <- rgb(0, 200, 190, alpha=50, maxColorValue=255)
  points(bigRes$emulatedx, bigRes$emulatedy+samples2, col=colTest, pch=16, type="p", cex=0.75 )
}
points(ourModel$xmodel, ourModel$training, ylim=range(0.0,6.0), xlab="x", ylab="y", pch=19)
title(main='Naiive Sample Cov Matrix')
lines(actual$x, actual$y, col="black", lwd=2, lty=2)
lines(results$emulatedx, results$emulatedy, col="red", lwd=2)
lines(results$emulatedx, results$emulatedy + confidence, col="red", lty=2)
lines(results$emulatedx, results$emulatedy - confidence, col="red", lty=2)

## legend(x=0.8, y=5, legend=c('model', 'emulator', 'confidence', 'samples badly'),
##        col=c('black', 'red', 'red', colTest),
##        lwd=2,
##        lty=c(2,1,2,1))


dev.off()
