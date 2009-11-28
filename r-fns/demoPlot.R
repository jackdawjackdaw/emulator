## a demo
source("EmuRbind.R") ## provides the basic interface to the c fns via libRbind
source("testRbind.R") ## provides some model test functions for plotting etc

# npts
npts <- 8


# this function in testRbind.R creates some model data
# with the yM wobbly sin(exp) fn over the range [0,1]
# the first arg sets the number of points to distribute in that range
# the second arg sets the if we should use lhs or not (we do here)
ourModel <- demoModel(npts, lhs=1)

# now we set the default values for the covariance function
# we'll need to know this later when we want to rebuild the
# C matrix
setDefaultOps()

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

## localCovFn
