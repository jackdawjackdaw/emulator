## some functions for testing the estimation parts of the emulator
##
## ccs, cec24@phy.duke.edu, aug-2011
##
## aims:
## 1: plot lhood over a plane while keeping other hyperparams fixed (lhoodcompare)
## 2: compare max lhood for various regression and cov fn models (not implemented yet)
## 3: do a k-sets withold test to build up a distribution of deviates (holdbacktest)
##
## inputs:
## modelData, trainingData.
## these will be supplied in model.in{training=trainingData, xmodel=modelData}
##
## these fns are for a single observable, we're not using multidim directly


##
## rangeMatrix: 2x2, { { minA, maxA}, {minB, maxB} }
## thetaVec, vector of the "best" thetas for this function
## we'll run thetaA and thetaB over the supplied ranges while holding the other values
## fixed
lhoodOverRange <- function(thetaA, thetaB, rangeMatrix, thetaVec, nthetas, nEvalPts, model.in)
{
  lhood.mat <- matrix(0, nrow=nEvalPts, ncol=nEvalPts)
  nmodelPoints <- dim(model.in$xmodel)[1]
  nparams <- dim(model.in$xmodel)[2]

  cat("nparams: ", nparams, "\n")
  cat("nmodelPoints: ", nmodelPoints, "\n")

  stepA <- abs((rangeMatrix[1,1] - rangeMatrix[1,2]) / nEvalPts)
  stepB <- abs((rangeMatrix[2,1] - rangeMatrix[2,2]) / nEvalPts)
  
  pointList <- matrix(0, nrow=nEvalPts**2, ncol=nthetas)
  for(i in 1:nEvalPts){
    for(j in 1:nEvalPts){
      subValA <- rangeMatrix[1,1] + i * stepA
      subValB <- rangeMatrix[2,1] + j * stepB
      subVec <- thetaVec
      subVec[thetaA] <- subValA
      subVec[thetaB] <- subValB
      pointList[j+nEvalPts*(i-1),] = subVec
    }
  }

  lhoodRes <- callEvalLhoodList(model.in, pointList, nEvalPts**2, nmodelPoints, nparams=nparams, nthetas=nthetas)

  lhoodRes
}

##
## vary 2 hyperparams and plot the lhood while keeping the remaining hyperparams fixed
## plot this as a set of images
lhoodCompare <- function(estim.result=estimResult, obsIndex=1, ngridPts=32, titleAdditional=NULL){
  namesfirst <- c("scale", "nugget")
  nthetas <- length(estim.result$thetas[obsIndex,])
  par(mfrow=c(nthetas,nthetas), mar=c(2,2,0,0), oma=c(2,2,2,2))

  # compute the lhood of the initial theta set
  model <- list(xmodel=estim.result$des.scaled, training=estim.result$train.scaled)
  lhoodMain <- callEvalLhoodList(model, estim.result$thetas[obsIndex,], 1,
                                 nmodelPoints=dim(model$xmodel)[1], nparams=dim(model$xmodel)[2], nthetas=nthetas)

  print(lhoodMain)
  
  
  # all this to make a grid
  for(i in 1:nthetas){
    for(j in 1:nthetas){
      # don't plot in the upper right corner
      if(j > i){
        plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
      } else if(j == i){
        # on the diagonals plot the names of the variables
        plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
        
        if(i <= 2){
          buffer <- namesfirst[i]
        } else {
          buffer <- paste("theta:", (i-2))
        }
        text(x=0.5, y=0.5, buffer, cex=1.5)
      } else {
        # actually plot a slice in the lhood plane
        plotLhoodSlice(estim.result, j, i, obsIndex, ngridPts)
      }
    }
  }

  par(mfrow=c(1,1), mar=c(2,2,2,2))
  buffer<-paste(titleAdditional, "L = ", round(lhoodMain$lhood,digits=2))
  title(main=buffer, line=-1)
}


plotLhoodSlice <- function(estim.result=estimResult, dimA=5, dimB=6, obsIndex=1, ngridPts=32, nconts=10)
{

  ranges <- matrix(0, nrow=2, ncol=2)
  if(dimA > 2){
    ranges[1,] <- c(-10,10)
  } else if (dimA == 2){
    ranges[1,] <- c(0.00001, 0.5)
  } else {
    ranges[1,] <- c(0.00001, 5)
  }

  if(dimB > 2){
    ranges[2,] <- c(-10,10)
  } else if (dimB == 2){
    ranges[2,] <- c(0.00001, 0.5)
  } else {
    ranges[2,] <- c(0.00001, 5)
  }

  # used the estimated thetas for the first bin
  fixedValVec <- estim.result$thetas[obsIndex,]
  print(fixedValVec)
  cat("ranges A: " , ranges[1,], "\n")
  cat("ranges B: " , ranges[2,], "\n")
  
  model <- list(xmodel=estim.result$des.scaled, training=estim.result$train.scaled[,obsIndex])
  
  res <- lhoodOverRange(dimA,dimB ,ranges, fixedValVec, nthetas=length(fixedValVec), nEvalPts=ngridPts, model)

  
  xrange <- seq(from=ranges[1,1], to=ranges[1,2], length.out=ngridPts)
  yrange <- seq(from=ranges[2,1], to=ranges[2,2], length.out=ngridPts)

  lhoodMatrix <- t(matrix(res$lhood, ncol=ngridPts, nrow=ngridPts))

  bufferX <- paste("theta :", dimA)
  bufferY <- paste("theta :", dimB)

  #levels=c(-300, -200,-150, -100, -75, -65,  -50, -40, -32, -25, -15, -10)
  #image(xrange, yrange, lhoodMatrix, xlab=bufferX, ylab=bufferY, col=topo.colors(length(levels)-1), breaks=levels)
  image(xrange, yrange, lhoodMatrix, xlab=bufferX, ylab=bufferY)
  contour(xrange, yrange, lhoodMatrix, nlevels=nconts, add=TRUE, labcex=0.8)
  points(fixedValVec[dimA], fixedValVec[dimB], pch=3, cex=1.5, lwd=2)

}

# partition the nmodelpoints into nSets subsets, estimate the model 
# for the complement of each subset and then compute the deviation on the
# witheld points
HoldBackTest <- function(des.scaled, train.scaled, nSets=7, fixNuggetOption=NULL){
  nmodelPts.tot <- dim(des.scaled)[1]
  nparams <- dim(des.scaled)[2]
  nthetas <- nparams + 2 ## currently
  # the number of points in one of the withold sets

  subSetLength <- floor( nmodelPts.tot / nSets)

  if(nmodelPts.tot %% nSets != 0){
    nSets <- nSets + 1 # add an additional set to hold the remainder
  }
  
  model.list <- vector("list", nSets)
  partial.des <- vector("list", nSets)
  partial.train <- vector("list", nSets)
  
  test.pts <- vector("list", nSets)
  true.values <- vector("list", nSets)

  subList <- vector("list", nSets)

  thetas.full <- matrix(0, ncol=nthetas, nrow=nSets)
  lhoods <- rep(NA, nSets)
  
  # loop over all the sets creating the partial designs, we make sure that the
  # final set contains the remainder
  for(i in 1:nSets){
    if(i < nSets){
      subList[[i]] <- seq(from=((i-1)*subSetLength+1), to=(i*subSetLength))
    } else {
      subList[[i]] <- seq(from=((i-1)*subSetLength+1), nmodelPts.tot)
    }
    # these will be the design points for each partial run
    partial.des[[i]] <- des.scaled[-(subList[[i]]), ]
    # the training values for each partial run
    partial.train[[i]] <- train.scaled[-(subList[[i]])]
    # the test locations
    test.pts[[i]] <- des.scaled[subList[[i]], ]
    # the values at the test locations
    true.values[[i]] <- train.scaled[subList[[i]]]

    model.list[[i]] <- list(xmodel=partial.des[[i]], training=partial.train[[i]])
  }

  # estimate the hyper params for each model 
  # we support the fixedNugget option here.
  #
  for(i in 1:nSets){
    npts <- length(model.list[[i]]$training)
    thetas.full[i,] <- callEstimate(model.list[[i]], nmodelpts=npts, nparams=nparams,
                               nthetas=nthetas, fixedNugget=fixNuggetOption)
  ## it's also fun to compute the lhood for each set
    temp <- callEvalLhoodList(model.list[[i]], thetas.full[i,], nevalPoints=1, nmodelPoints=npts, 
                                   nparams=nparams, nthetas=nthetas)
    lhoods[i] <- temp$lhood
  }

  # hold the results and the deviates
  emu.res <- vector("list", nSets)
  devs <- c()
  
  # now we want to compute the deviates for each set
  for(i in 1:nSets){
    nEmuPts <- length(true.values[[i]])
    npts <- length(model.list[[i]]$training)
    # start by emulating each partial model at its test.pts locations
    emu.res[[i]] <- callEmulateAtList(model.list[[i]], thetas.full[i,], test.pts[[i]], nemupts=nEmuPts,
                              nmodelpoints=npts, nparams=nparams, nthetas=nthetas)

    ## we quantify the deviations as the squared diff from the real value scaled by 
    ## the predicted variance at this locn
    devs <- c(devs, (emu.res[[i]]$mean - true.values[[i]])**2 / (emu.res[[i]]$var))
  }

  ## we return the vector of deviates, the test points where we compute the deviates
  ## the emulated means and vars at these test points, the true values at the points
  ## the lhoods of each of the models and the thetas specifying each model
  ## and also the models themselves
  final.res <- list(deviates=devs, test.pts=test.pts, emu.res=emu.res,
                    true.values=true.values, lhoods=lhoods, thetas=thetas.full,
                    models=model.list)
  
}

# use this to test holdbacktest
testWT <- function(){
  options(error=utils::recover) 
  source("emuYields.R")
  res <- HoldBackTest(estimResult$des.scaled, estimResult$train.scaled[,1], nSets=20)
  invisible(res)
}



