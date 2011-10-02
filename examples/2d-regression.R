## regression in 2d parameter space
## ccs, cec24@phy.duke.edu
##
## an example of using the R interface to make predictions of a 2d function
##
## we follow the same basic procedure as in simple-regression.R:
## - generate a set of samples from a model
## - estimate appropriate thetas for the model
## - make predictions on a grid of untried points
## - plot our predictions and compare against the full model
##
## the only tricky part with going to higher dimensions (apart from visualising the result)
## is making sure that libRBIND understands what you've sent it, in the R terminal
## you should see lines like: nx = 2, ny = 32
## where nx is the dimensionality of your design space and ny is the number of training points
## if you're seeing nx = 32 and ny =2 then either your xmodel is transposed, it should have
## nmodelpoints rows and nparams columns or your definitions of nparams and nmodelpoints are
## broken somewhere
## 

## load EmuRbind
## 
source("~/local/include/libRbind/EmuRbind.R") 

## load libRBIND into R
arch <- system("uname -s", intern=TRUE)
if(is.loaded("callEstimate") == FALSE){
  libNAME <- "~/local/lib/libRBIND"
  if(arch == "Linux"){
    libNAME <- paste(libNAME, ".so", sep="")
  } else if(arch =="Darwin"){
    libNAME <- paste(libNAME, ".dylib", sep="")
  } else {
    buffer <- paste("error: uname -s gives ", arch , " not supported", sep="")
    stop(buffer)
  }
  dyn.load(libNAME)
}

library(mvtnorm)

##
## some wobbly gaussians in 2d
yM <- function(vec, delta=0.0005)
{
  sig1<-  matrix(0, nrow=2, ncol=2)
  sig1[1,1] <- 0.4
  sig1[1,2] <- 0.2
  sig1[2,2] <- 0.8
  sig1[2,1] <- 0.2
  f <- dmvnorm(vec, mean=c(0.5, 1.3), sigma=0.1*sig1) + rnorm(1,mean=0, sd=delta) +
    dmvnorm(vec, mean=c(1.2, 1.8), sigma=0.1*diag(2)) + 
      dmvnorm(vec, mean=c(1.58, 0.9), sigma=0.2*sig1) +
   dmvnorm(vec, mean=c(0.7, 0.8), sigma=0.2*diag(2))

  f
}

## construct the model data for emulation using a latin hypercube sampling method to generate the design
## @return a list containing the model design and the evaluations of the model over the design
make.model.2d <- function(m, rangeMin=0.0, rangeMax=2.0, modelFunc=yM){
  nparams <- 2
  retval <- require("lhs", quietly=TRUE)

  ## test for the library
  if(retval == FALSE){
    print("installing the lhs library")
    install.packages(c("lhs"))
  }
  ## we use the R library LHS to construct a latin hypercube design in 2d
  design <- (maximinLHS(m, nparams))*(rangeMax-rangeMin) + rangeMin

  ymodel <- rep(NA, m)
  for(i in 1:m)
    ymodel[i] <- modelFunc(design[i,])
  
  ## the model canonically contains an xmodel and a training set list, other entries
  ## are ignored
  model  <- list(xmodel = design, training = ymodel)
  invisible(model)
}


## estimate the hyperparameters for the supplied model
## @return a list containing the model, the estimated thetas and all the parameters needed to
## generate predictions from the model at some new untried input
estimate.model <- function(model){
  nmodelpts <- dim(model$xmodel)[1]   ## number of design points (number of rows in xmodel)
  nparams <- 2 ## number of dimensions in the u space this is now 2!
  ## for the default power-exponential covariance function this nthetas = nparams + 2
  ## callEstimate will complain if this is set to an incorrect value.
  nthetas <- nparams + 2
  cov.fn <- 1
  reg.order <- 1

  ## callEstimate uses libRBIND to try and generate the best set of length scales for the supplied model
  thetas.est <- callEstimate(model, nmodelpts,
                             nparams=nparams, nthetas=nthetas,
                             cov.fn=cov.fn, reg.order=reg.order)

  ## return a list with the appropriate thetas, the model data used to generate them
  estim.result <- list(thetas=thetas.est, model=model, cov.fn=cov.fn, reg.order=reg.order,
                       nparams=nparams, nmodelpts=nmodelpts, nthetas=nthetas)
  invisible(estim.result)
}


## generate predictions for the model using the results of estimate.model
## @param estim.result list returned from estimate.model
## @return a list the points we evaluated the model at and the predicted mean and variance at each of these
predict.model <- function(estim.result){
  nemupts <- 64

  ## generate a uniform square grid of points across our 2d space
  xmax <- max(estim.result$model$xmodel[,1])
  xmin <- min(estim.result$model$xmodel[,1])
  x <- seq(from=xmin, to=xmax, length.out=nemupts)
  pointList <- as.matrix(expand.grid(x, x))

  
  ## callEmulateAtList uses libRBIND to generate predictions for the supplied model at
  ## the set of pointList locations in the parameter space.
  ## 
  ## we specify the model with the results from estimate.model
  ##
  ## note that we supply nemupts**2 since we're sampling on a 64 x 64 grid
  emu.result <- callEmulateAtList(estim.result$model, estim.result$thetas,
                                  pointList, nemupts**2, nmodelpoints=estim.result$nmodelpts,
                                  nparams=estim.result$nparams, nthetas=estim.result$nthetas,
                                  reg.order=estim.result$reg.order,
                                  cov.fn=estim.result$cov.fn)
  invisible(emu.result)
}


## plot the results
##
## we show the mean and variance separately and then compute a
## variance weighted distance of the true function from the
## emulator mean
plot.result.2d <- function(emu.result, model, modelfn){
  npts <-  64
  par(mfrow=c(2,2))

  xmax <- max(model$xmodel[,1])
  xmin <- min(model$xmodel[,1])
  x <- seq(from=xmin, to=xmax, length.out=npts)

  mean.mat <- matrix(emu.result$mean, npts, npts)
  var.mat <- matrix(emu.result$var, npts, npts)

  ## first plot the mean of the emulator
  image(x,x,mean.mat)
  contour(x,x,mean.mat, add=TRUE)
  title(main="emulated mean")

  ## the varianace
  image(x,x,var.mat)
  contour(x,x,var.mat, add=TRUE)
  title(main="emulated variance")

  ## evaluate the original model over a uniform grid
  ## this is actually a lot slower than evaluating the emulator over
  ## the same grid...
  grid.orig <- expand.grid(x,x)
  grid.orig.eval <- apply(grid.orig, 1, modelfn)
  image(x,x,matrix(grid.orig.eval,npts,npts))
  contour(x, x, matrix(grid.orig.eval, npts, npts), add=TRUE)
  points(model$xmodel[,1], model$xmodel[,2], pch=2, col="green")
  title(main="original function")

  ## the deviation
  dev <- sqrt((mean.mat - grid.orig.eval)**2 / (var.mat))
  image(x,x,dev, col=cm.colors(16))
  contour(x, x, dev, add=TRUE)
  title(main="deviation from original")
  

}

## sample our original function with a lhs design
## changing this argument is interesting, < 32 gives rather poor (or variable) results
## the computational process starts to get quite slow for ~ 100 points
##
model <- make.model.2d(256)

## estimate optimal hyper-parameters for this model
estimResult <- estimate.model(model)
## make predictions over a uniform 2d grid
emuResult <- predict.model(estimResult)
## plot the predicted mean and var along with the original model
##
## it's interesting to examine the deviation structure, our original functions have
## some non trivial x~y variation and this seems to come out for nmodelpoints ~ 100
## as small ellipses aligned along x~y which are not well treated
##
## pushing up to ~ 200 fixes this (huge oversample)
plot.result.2d(emuResult, model, yM)
