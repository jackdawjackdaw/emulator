## ccs,
## make a simple emulator and plot samples of it
## check that these samples match up with our confidence intervals
##
## the samples seem to now massively underestimate the confidence intervals
## my guess is that this is a result of the regression process which adds a
## term to the covariance

arch <- system("uname -s", intern=TRUE)
if(arch == "Linux"){ 
  dyn.load("~/local/lib/libRBIND.so") # load the emulator
} else if (arch == "Darwin"){
  dyn.load("~/local/lib/libRBIND.dylib") # load the emulator
} else {
  buffer <- paste("error: uname -s gives ", arch , " not supported", sep="")
  stop( buffer )
}
source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings
if(is.loaded("callEstimate") == FALSE){
  stop("loading of libRBIND didn't work")
}

## for making some demo plots
yM <- function(z) {5*exp(-3*z)*(sin(z*10)) + 2}

## the data set for this
demoModel <- function(m, lhs=0, rangeMin=0.0, rangeMax=1.0, modelFunc=yM){
  if(lhs == 1){
    params <- (maximinLHS(m, 1))*(rangeMax-rangeMin) + rangeMin
  } else {
    params <- matrix(0,m,1)
    x<-  seq(from=rangeMin,to=rangeMax,length.out=m)
    for(i in 1:m)
      params[i,] <- x[i]

  }
  ymodel = rep(NA, m)
  for(i in 1:m)
    ymodel[i] <- modelFunc(params[i,])
  model <- data.frame(xmodel = params, training = ymodel)
  model
}


emulateSimpleModel <- function(){
  nemupts <- 128
  
  #input.data <- as.matrix(read.table("../input/model-cut.dat"))
  #design <- input.data[,1]
  #training <- input.data[,2]
  #nmodelpts <- length(design)
  #model <- list(xmodel=scale(design), training=scale(training))
  nmodelpts <- 10
  model <- demoModel(nmodelpts, lhs=1)

  estim.result <- callEstimate(model, nmodelpts, nparams=1, nthetas=3)

  pointList <- seq(from=min(model$xmodel), to=max(model$xmodel), length.out=nemupts)

  emu.result <- callEmulateAtList(model, estim.result, pointList, nemupts, nmodelpts, nparams=1, nthetas=3)

  finalVal <- list(thetas=estim.result, emu.result=emu.result, model=model)

  invisible(finalVal)
}

## generate samples, we're not using the regression model right now
sampleModel <- function(emu.result, nsamples=128, nsamplePts=128){
  des.orig <- emu.result$model$xmodel
  nmodelpts <-length(des.orig)
  sampleLocs <- seq(from=min(des.orig), to=max(des.orig), length.out=nsamplePts)

  # this is the mean vector for all the samples
  sampleMean <- callEmulateAtList(emu.result$model, emu.result$thetas, sampleLocs, nemupts=nsamplePts,
                                  nmodelpoints=length(des.orig), nparams=1, nthetas=3)

  # this will be our random vector
  sample.rand <- rnorm(nsamplePts)

  # start by evaluating the covariance of each of the new locations
  cov.matrix <- matrix(0, nrow=nsamplePts, ncol=nsamplePts)
  for(i in 1:nsamplePts){
    for(j in 1:nsamplePts){
      cov.matrix[i,j] <- cov.fn(sampleLocs[i], sampleLocs[j], emu.result$thetas)
      if(i == j){
        cov.matrix[i,j] <- cov.matrix[i,j] + emu.result$thetas[2] # add in the nugget
      }
    }
  }

  # now we eval the original cov matrix (ish)
  cov.matrix.orig <- matrix(0, nrow=nmodelpts, ncol=nmodelpts)
  for(i in 1:nmodelpts){
    for(j in 1:nmodelpts){
      cov.matrix.orig[i,j] <- cov.fn(des.orig[i], des.orig[j], emu.result$thetas)
    }
  }

  # epsilon is here to add some numerical stability
  epsilon <- 0.0001
  diag(cov.matrix.orig) <- diag(cov.matrix.orig) + epsilon
  
  cov.matrix.orig.inv <- chol2inv(chol(cov.matrix.orig))

  # now we make the marginal K matrix (not a vector)
  # this is the matrix of the covariance of each sample point
  # against each of the original model points
  k.matrix <- matrix(0, nrow=nsamplePts, ncol=nmodelpts)
  for(i in 1:nsamplePts){
    for(j in 1:nmodelpts){
      k.matrix[i,j] <- cov.fn(sampleLocs[i], des.orig[j], emu.result$thetas)
    }
  }

  # our final covariance structure is the covariance of the sample locations
  # reduced by the training data

  sigma <- cov.matrix - (k.matrix %*% cov.matrix.orig.inv) %*% t(k.matrix)
 
  sample <- sampleMean$mean + sigma %*% sample.rand
  
  result <- list(x=sampleLocs, y=sample)
}

cov.fn <- function(x, y, thetas){
  lscale <- exp(thetas[3])
  ans <- thetas[1] * exp(-(x-y)^2/lscale)
}

nsamples <- 256

emu.result <- emulateSimpleModel()
samp.list <- vector("list", nsamples)

plot(emu.result$emu.result$des, emu.result$emu.result$mean, type="l", col="red", ylim=c(-5,5))
points(emu.result$model$xmodel, emu.result$model$training, pch=3, lwd=2)


for(i in 1:nsamples){
  samp.list[[i]] <- sampleModel(emu.result, nsamples=1, nsamplePts=128)
  colTest <- rgb(0, 102, 153, alpha=35, maxColorValue=255)
  lines(samp.list[[i]]$x, samp.list[[i]]$y, col=colTest)
}

confUp <- emu.result$emu.result$mean + 1.64*sqrt(emu.result$emu.result$var)
confDown <- emu.result$emu.result$mean - 1.64*sqrt(emu.result$emu.result$var)
lines(emu.result$emu.result$des, confUp, col="red", lwd=2, lty=2)
lines(emu.result$emu.result$des, confDown, col="red", lwd=2, lty=2)
