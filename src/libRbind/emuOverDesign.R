## ccs, cec24@phy.duke.edu
## functions for emulating scalar model output (potentially many independent scalars)
## over a design.
##
## Aug-2011
## \todo add a gridImplausSweep
## \todo add a implausObsOverDesign fn


## doCombEstimation -> computes optimal hyperparameters and produces a result list
##
## inputs (global):
## designData, a matrix rows of which are the design location of a
## single run, there should be nruns rows and nparams columns
##
## modelData, a matrix rows of which are the model output for each of the scalar observables.
## there should be nruns rows and nobservables parameters
##
## the design and model output are centered and scaled, means subtracted and then divided by their
## sample variance. This scaled data set is used to train a set of independent emulators through
## the multidim call. 
## 
##
## outputs:
## a list with fields, thetas, des.scaled, train.scaled
## thetas  -> a matrix of optimal hyperparams, there are nobservable rows each with nparams columns
## des.scaled -> the scaled and centered design
## train.scaled -> the scaled and centered training points (code output)
## 
## 
doCombEstimation <- function(){
# now we want to scale the design onto a unit cube, otherwise the estimation etc is not going to work well
  scaledDesign <- scale(designData)
# do we want to scale the output too?
  scaledModelData <- scale(modelData) # well this does make things nicer
  combinedThetas <- multidim(scaledDesign, nmodelpts=nruns, training=scaledModelData, nydims=nbins)
  result <- list(thetas=combinedThetas, des.scaled=scaledDesign, train.scaled=scaledModelData)
  invisible(result)
}

##
## emulates the given observable over the ranges of design vectors A and B
## 
## we keep the other observables fixed at the values given in the  "fixedValVec"
## the thetas should just be the vector of thetas for this observable
## nEmuPts is the number of points per design variable (the sqrt of the total)
## 
## fullDesign is there as a kludge
emulateObsOverDesign <- function(observable, thetas, fullDesign,
                                 designA, designB, fixedValVec,nEmuPts){
  rangeA <- c(min(designA), max(designA))
  rangeB <- c(min(designB), max(designB))
  nparams <- 2 + length(fixedValVec)
  nmodelpoints <- length(designA)
  
  stepSizeA <- (rangeA[2] - rangeA[1]) / (nEmuPts)
  stepSizeB <- (rangeA[2] - rangeA[1]) / (nEmuPts)


  
  mean <- matrix(0, nrow=nEmuPts, ncol=nEmuPts)
  var <- matrix(0, nrow=nEmuPts, ncol=nEmuPts)
  pointList <- matrix(0, nrow=nEmuPts**2, ncol=nparams)

  browser()
  
  for(i in 1:nEmuPts){
    for(j in 1:nEmuPts){
      pointList[j+nEmuPts*(i-1),] = rbind(rangeA[1]+i*stepSizeA, rangeB[1]+j*stepSizeB, fixedValVec)
    }
  }

 # browser()
  
  # now we can emulate all these points at once
  emuRes <- callEmulateAtList(list(xmodel=fullDesign, training=observable),
                              thetas, pointList,
                              nemupts = dim(pointList)[1],
                              nmodelpoints = nmodelpoints,
                              nparams = nparams, # at least A and B, right?
                              nthetas = nparams + 2) # for the power exp model
  emuRes
}

## 
## emulates the observable indicated over the ranges of design variable A and B with any remaining variables
## being held at the values in fixedVal.
## 
## All the model data comes in through estim.result
##
## output: data is plotted on the current device
##
## inputs:
## obsIndex -> index of the observable we care about in fullReslut
## dAIndex, dBIndex -> index of the two design variables we'll span
## fixedVal -> vector of the remaining (nparams - 2) values which are held fixed while we vary the A and B
## xlabel, ylabel -> human-readable names of dA and dB. (goes on the plot)
## plotDes -> if true plots the design points projected onto A,B
## unscale -> if true undoes scaling and centering (plot values should be back in reasonable range)
## estim.result -> list of thetas, des.scaled, train.scaled created by doCombEstimation
## 
plotEmuOverDesign <- function(obsIndex, dAIndex, dBIndex, fixedVal,
                              xlabel, ylabel, titleIn, plotDes=FALSE, unscale=TRUE, estim.result){
  npts <- 32
  print(estim.result$thetas[obsIndex,])

  #browser()
  designA <- estim.result$des.scaled[,dAIndex]
  designB <- estim.result$des.scaled[,dBIndex]
  
  emu.result <- emulateObsOverDesign(estim.result$train.scaled[,obsIndex], estim.result$thetas[obsIndex,],
                                     estim.result$des.scaled, designA,
                                     designB, fixedVal, npts)

  
  mean <- matrix(emu.result$mean, nrow=npts, ncol=npts)
  var <- matrix(emu.result$var, nrow=npts, ncol=npts)

  range1 <- seq(min(designA), max(designA), length=npts)
  range2 <- seq(min(designB), max(designB), length=npts)
  
  if(unscale==TRUE){ # we'll return to the original scale for plotting
    trainCenter <- attr(estim.result$train.scaled, "scaled:center")[obsIndex]
    trainScale <- attr(estim.result$train.scaled, "scaled:scale")[obsIndex]

    mean <- mean*trainScale + trainCenter
    # V(a*x) = a**2 V(x)
    # man i suck
    var <- var*(trainScale**2)

    desACenter <- attr(estim.result$des.scaled, "scaled:center")[dAIndex]
    desBCenter <- attr(estim.result$des.scaled, "scaled:center")[dBIndex]

    desAScale <- attr(estim.result$des.scaled, "scaled:scale")[dAIndex]
    desBScale <- attr(estim.result$des.scaled, "scaled:scale")[dBIndex]

    range1 <- range1 * desAScale + desACenter
    range2 <- range2 * desBScale + desBCenter
        
  }

  
  image(range1, range2, t(mean), axes=FALSE, col=heat.colors(16), xlab=xlabel, ylab=ylabel)
  contour(range1, range2, t(mean), nlevels=10, col="black", add=TRUE, cex.lab=0.5, labcex=0.8)
  title(main=titleIn)
  if(plotDes==TRUE){
    points(designA, designB, pch=3)
  }
  axis(1, cex.axis=2)
  axis(2, cex.axis=2)
}

## midPoint <- mean(estimResult$des.scaled[,3])
## plotEmuOverDesign (estimResult$train.scaled[,1], estimResult$thetas[1,], estimResult$des.scaled[,2], estimResult$des.scaled[,3], c(midPoint), "", "", "")
                     

##
## repeatedly calls plotEmuOverDesign creating a set of 9 plots of obsIndex which span the stepping dimension stepDim
##
## fixedVals -> a vector of additional fixed values (empty unless nparams > 3)
stepPlotDimension <- function(estim.result=estimResult, obsIndex, plotDimA, plotDimB, stepDim, fixedVals=NULL, nsteps=9){

  minVal <- min(estim.result$des.scaled[,stepDim])
  maxVal <- max(estim.result$des.scaled[,stepDim])
  
  stepSize <- (maxVal - minVal) / nsteps

  par(mfrow=c(3,3))

  for(i in 0:(nsteps-1)){
    fixV <- minVal + stepSize * i

    # now if we're unscaling we need to unscale the fixed value also
    fixScale <- attr(estim.result$des.scaled, "scaled:scale")[stepDim]
    fixCenter <- attr(estim.result$des.scaled, "scaled:center")[stepDim]
    fixVUnscaled <- fixV * fixScale + fixCenter
    
    buffer <- paste(desNames[stepDim], " fixed at: ", round(fixVUnscaled,2), sep="")
    if(i == 0){
      desPlot <- TRUE
    } else {
      desPlot <- FALSE
    }
    plotEmuOverDesign(obsIndex, plotDimA ,plotDimB, c(fixV, fixedVals),
                      xlabel=desNames[plotDimA], ylabel=desNames[plotDimB], titleIn=buffer,
                      desPlot, unscale=TRUE, estim.result)
  }
}

##
## run the emulator over variables A,B,C and output the resulting "huge" grid
##
## by default we unscale the emulator results
##
## 
gridEmulatorSweep <- function(estim.result=estimResult, obsIndex, dimA, dimB, dimC, fixedVals=NULL, nGridPts = 32,
                              unscale=TRUE, fname="grid-sweep.dat"){

  ndim <- 3
  nEmuPts <- nGridPts**ndim

  des <- estim.result$des.scaled
  obs <- estim.result$train.scaled[,obsIndex]

  nParams <- dim(des)[2]
  nModelPoints <- dim(des)[1]

  pointList <- matrix(0, nrow=nEmuPts, ncol=nParams)
  ranges <- matrix(0, nrow=ndim, ncol=2)

  dimList <- c(dimA, dimB, dimC)
  # read the names from the design
  dimnames <- attr(estim.result$des.scaled, "dimnames")[[2]]
  axisNames <- rep(NA, ndim)
  for(i in 1:ndim){
    axisNames[i] <- dimnames[dimList[i]]
  }
  # contains the explicit axes, all the points will be some combination of these ndim column vectors
  # we can use expand.grid to make this combination for us
  axes <- matrix(0, nrow=nGridPts, ncol=ndim)

  for(i in 1:ndim){
 
    ranges[i,] <- c(min(des[,dimList[i]]), max(des[,dimList[i]]))
    
    axes[,i] <- seq(from=ranges[i,1], to=ranges[i,2], length.out=nGridPts)
  }

  cat("nparams: ", nParams, " \n")
  cat("nmodelpts: ", nModelPoints, " \n")

  # can't figure out how to do this in a dimensionally agnostic way here..?
  pointGrid <- expand.grid(A=axes[,1], B=axes[,2], C=axes[,3])

  if(is.null(fixedVals) != TRUE){
    fixedValsSeq <- matrix(0, nrow=nEmuPts, ncol=dim(fixedVals))
    for(i in 1:nEmuPts){
      fixedValSeq[i,] <- fixedVals
    }
    pointList <- as.matrix(cbind(pointGrid, fixedValSeq))
  } 
  
  pointList <- as.matrix(pointGrid)
  
  # now we emulate all these points (this will probbaly be slow)
  emuRes <- callEmulateAtList( list(xmodel=des, training=obs),
                              estim.result$thetas[obsIndex,], pointList,
                              nemupts = nEmuPts, nmodelpoints = nModelPoints,
                              nparams = nParams, nthetas = nParams + 2 )

  # now unscale the data
  if(unscale==TRUE){
    centerScales.Des <- matrix(0, nrow=ndim, ncol=2)
    for(i in 1:ndim){
      centerScales.Des[i,] <- c( attr(estim.result$des.scaled, "scaled:center")[dimList[i]],
                                  attr(estim.result$des.scaled, "scaled:scale")[dimList[i]] )

      # for each axis(column), we undo the scaling by multiplying by the removed sample variance (the scale)
      # and adding back in the mean (the center)
      pointList[,i] <- pointList[,i] * centerScales.Des[i, 2] + centerScales.Des[i,1] 
    }
    centerScales.Train <- c( attr(estim.result$train.scaled, "scaled:center")[obsIndex],
                                  attr(estim.result$train.scaled, "scaled:scale")[obsIndex] )
    
    mean <- emuRes$mean * centerScales.Train[2] + centerScales.Train[1]
    var <- emuRes$var * (centerScales.Train[2]**2)
  } else {
    mean <- emuRes$mean
    var <- emuRes$var
  }

  browser()
  
  finalTable <- cbind(pointList, mean, var)

  # we write the table out to a csv file for vis
  write.table(finalTable, file=fname, row.names=FALSE, col.names=c(axisNames, "mean", "variance"),
              sep=",", qmethod="double", dec=".")

  invisible(finalTable) # returns final table if you want it
  
}
