## a demo
source("EmuRbind.R") ## provides the basic interface to the c fns via libRbind
source("testRbind.R") ## provides some model test functions for plotting etc
source("plotCovReals.R") ## makes a covariance matrix
library("Matrix")



# npts
npts <- 7


# this function in testRbind.R creates some model data
# with the yM wobbly sin(exp) fn over the range [0,1]
# the first arg sets the number of points to distribute in that range
# the second arg sets the if we should use lhs or not (we do here)
ourModel <- demoModel(npts, lhs=0)

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
results <- callEmulate(ourModel, thetas, npts, rangemin=0.0, rangemax=1.5)

bigRes <- callEmulate(ourModel, thetas, npts, nemupts=400, rangemin=0.0, rangemax=1.1)

# this is going to be our analytic model for plotting against
# not transparent but the fn is def'd in testRbind
sequence <- seq(0.0, 1.0, length=100)
actual <- data.frame(x=sequence, y=yM(sequence))


# we can plot the results, and it works!
# plot a jpg
#jpeg("model-showing-samples.jpg", quality=100, bg="white", res=200, width=8, height=5, units="in")
#plotResultsTest(ourModel, results, actual, "our model")


##
## this is ok but really we should do this at quite a few points
## 

# now make the covMatrix
cM <- makeCMatrix(npts, thetas, ourModel$xmodel)
bigpts<- 400
#cMHuge <- makeCMatrix(bigpts, thetas, results$emulatedx)
cMHuge <- makeCMatrix(bigpts, thetas, bigRes$emulatedx)

# make the cholesky decomp,
# now have f2' . f2 = cM) up to numeric errors
f2 <- chol(cMHuge)

plot(actual$x, actual$y, type="n", xlim=range(0,1.0), ylim=range(0.0,6.0), xlab='x', ylab='y')

for(i in 1:20){
  z1 <- rnorm(bigpts)
# now we have some samples with the right correlation
  samples2 <- t(f2) %*% z1
  colTest <- rgb(190, 190, 190, alpha=30, maxColorValue=255)
  points(bigRes$emulatedx, bigRes$emulatedy+samples2, col=colTest, pch=16, type="p", cex=0.5 )

}

                                        # this works in 1d now
points(ourModel$xmodel, ourModel$training, ylim=range(0.0,6.0), xlab="x", ylab="y", pch=19)
title(main='Showing the emulator with samples drawn from the resulting distribution')
lines(actual$x, actual$y, col="green", lwd=2, lty=2)
lines(results$emulatedx, results$emulatedy, col="red", lwd=2)
confidence <- rep(NA, length(results$emulatedvar))
for(i in 1:length(results$emulatedvar))
  confidence[i] <- 0.5*sqrt(results$emulatedvar[i])*1.65585

lines(results$emulatedx, results$emulatedy + confidence, col="red", lty=2)
lines(results$emulatedx, results$emulatedy - confidence, col="red", lty=2)

## add a legend
legend(x=0.8, y=5, legend=c('model', 'emulator', 'confidence', 'samples'),
       col=c('green', 'red', 'red', colTest),
       lwd=2,
       lty=c(2,1,2,1))
grid()

dev.off()
