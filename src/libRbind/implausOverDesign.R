## ccs, cec24@phy.duke.edu
## functions for evaluating the implausibility of a set scalar model output
## over the experimental design
##
## Aug-2011
##
## uses functions from emuOverDesign


## we require the experimental data be in the form
## of a list with entries: obsValue; obsError

## 
## compute the implausibility over a pair of parameters A, B with the remaining parameters
## held at the vector fixedVal
## 
implausObsOverDesign <- function(obsIndex, dAIndex, dBIndex, fixedVal, estim.result, expData, npts=32){

  nEmuPts = npts**2
  obs <- estim.result$train.scaled[,obsIndex]
  thetas <- estim.result$thetas[obsIndex,]
  desA <- estim.result$des.scaled[,dAIndex]
  desB <- estim.result$des.scaled[,dAIndex]
  
  nparams <- dim(estim.result$des.scaled)[2]
  nModelPoints <- dim(estim.result$des.scaled)[1]

  # run the emulator over the 2d plane of desA, desB
  emu.result <- emulateObsOverDesign(obs, thetas, estim.result$des.scaled, desA, desB, fixedVal, npts)

  emu.mean <- matrix(emu.result$mean, nrow=npts, ncol=npts)
  emu.var <- matrix(emu.result$var, nrow=npts, ncol=npts)
  range1 <- seq(min(desA), max(desA), length=npts)
  range2 <- seq(min(desB), max(desB), length=npts)

  # do the unscaling
  trainCenterScale <- c(  attr(estim.result$train.scaled, "scaled:center")[obsIndex], 
                        attr(estim.result$train.scaled, "scaled:scale")[obsIndex])
  emu.mean <- emu.mean*trainCenterScale[2] + trainCenterScale[1]
  # the V(a*x) = a**2 * V(x)
  emu.var <- emu.var*(trainCenterScale[2]**2)

  
  desACenterScale <- c( attr(estim.result$des.scaled, "scaled:center")[dAIndex],
                       attr(estim.result$des.scaled, "scaled:scale")[dAIndex])
  desBCenterScale <- c( attr(estim.result$des.scaled, "scaled:center")[dBIndex],
                       attr(estim.result$des.scaled, "scaled:scale")[dBIndex])
  range1 <- range1 * desACenterScale[2] + desACenterScale[1]
  range2 <- range2 * desBCenterScale[2] + desBCenterScale[1]
  
  final.emu.data <- list(mean=emu.mean, var=emu.var, rA=range1, rB=range2, fixed=fixedVal)

  implaus <- computeImplaus(obsIndex, final.emu.data, expData)

  result <- list(feas=implaus, emu.data=final.emu.data)

  invisible(result)
}



# plot a 2d implausibility matrix with all nice things
# implausResult is the mega list obtained from implausObsOverDesign
plotImplausOverDesign <- function(obsIndex, dAIndex, dBIndex, fixedVal,
                        xlabel="", ylabel="", titleIn="", plotDes=FALSE, estim.result, exp.data){
  
  implaus.result <- implausObsOverDesign(obsIndex, dAIndex, dBIndex, fixedVal, estim.result, exp.data, npts=32)
    
  r1 <- implaus.result$emu.data$rA
  r2 <- implaus.result$emu.data$rB

  desA <- estim.result$des.scaled[,dAIndex]
  desB <- estim.result$des.scaled[,dBIndex]

  desACenterScale <- c(attr(estim.result$des.scaled, "scaled:center")[dAIndex],
                       attr(estim.result$des.scaled, "scaled:scale")[dAIndex])
  
  desBCenterScale <-c(attr(estim.result$des.scaled, "scaled:center")[dBIndex],
                      attr(estim.result$des.scaled, "scaled:scale")[dBIndex])

  desA <- desA * desACenterScale[2] + desACenterScale[1]
  desB <- desB * desBCenterScale[2] + desBCenterScale[1]

  image(r1, r2, t(implaus.result$feas), axes=FALSE, col=cm.colors(16), xlab=xlabel, ylab=ylabel)
  contour(r1, r2, t(implaus.result$feas), nlevels=10, col="black", add=TRUE, cex.lab=0.5, labcex=0.8)
  title(main=titleIn)

  if(plotDes==TRUE){
    points(desA, desB, pch=3)
  }

  axis(1, cex.axis=2)
  axis(2, cex.axis=2)

}

##
## repeatedly calls plotImplaus creating a set of 9 plots of obsIndex which span the stepping dimension stepDim
##
## fixedVals -> a vector of additional fixed values (empty unless nparams > 3)
stepPlotDimensionImplaus <- function(obsIndex, plotDimA, plotDimB, stepDim, fixedVals=NULL, nsteps=9,
                                     estim.result, exp.data){

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
    plotImplausOverDesign(obsIndex, plotDimA ,plotDimB, c(fixV, fixedVals),
                      xlabel=desNames[plotDimA], ylabel=desNames[plotDimB], titleIn=buffer,
                      desPlot,  estim.result, exp.data)
  }
}


# I = (EmuMean  - ExpValue)**2 / (Var_Model + Var_Emulator + Var_Data)
# for facundos data we have no var_model
# 
# here we're supposing that expData$obsError is in the form of a standard error
# should support including errors from the model too
computeImplaus <- function(obsIndex, emu.data, expData){
  num <- (emu.data$mean - expData$obsValue[obsIndex])**2

  if(is.null(expData$errModel) == TRUE && is.null(expData$errModel.sys)== TRUE){ 
    denom <- (emu.data$var + (expData$obsError[obsIndex])**2)
  } else {
    denom <- (emu.data$var + (expData$obsError[obsIndex])**2 + (expData$errModel[obsIndex]**2))
  }

  # scalar matrix divide i hope?
  I <- num / denom

  I
}


##
## compute the implaus  over variables A,B,C and output the resulting "huge" grid
##
## by default we unscale the emulator results
## almost verbatim the fn gridEmulatorSweep from emuOverDesign
## 
gridImplausSweep <- function(estim.result=estimResult, exp.data, obsIndex, dimA, dimB, dimC, fixedVals=NULL, nGridPts = 32,
                              unscale=TRUE, fname="grid-sweep-implaus.csv"){

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

  pointGrid <- expand.grid(A=axes[,1], B=axes[,2], C=axes[,3])

  if(is.null(fixedVals) != TRUE){
    fixedValsSeq <- matrix(0, nrow=nEmuPts, ncol=dim(fixedVals))
    for(i in 1:nEmuPts){
      fixedValSeq[i,] <- fixedVals
    }
    pointList <- as.matrix(cbind(pointGrid, fixedValSeq))
  } 
  
  pointList <- as.matrix(pointGrid)
  
  emuRes <- callEmulateAtList( list(xmodel=des, training=obs),
                              estim.result$thetas[obsIndex,], pointList,
                              nemupts = nEmuPts, nmodelpoints = nModelPoints,
                              nparams = nParams, nthetas = nParams + 2 )

  # now unscale the data
  # actually because paraview is retarded, we're not going to do this
  ## centerScales.Des <- matrix(0, nrow=ndim, ncol=2)
  ## for(i in 1:ndim){
  ##   centerScales.Des[i,] <- c( attr(estim.result$des.scaled, "scaled:center")[dimList[i]],
  ##                             attr(estim.result$des.scaled, "scaled:scale")[dimList[i]] )

  ##                                       # for each axis(column), we undo the scaling by multiplying by the removed sample variance (the scale)
  ##                                       # and adding back in the mean (the center)
  ##   pointList[,i] <- pointList[,i] * centerScales.Des[i, 2] + centerScales.Des[i,1] 
  ## }

  centerScales.Train <- c( attr(estim.result$train.scaled, "scaled:center")[obsIndex],
                          attr(estim.result$train.scaled, "scaled:scale")[obsIndex] )

  emu.mean <- emuRes$mean * centerScales.Train[2] + centerScales.Train[1]
  emu.var <- emuRes$var * (centerScales.Train[2]**2)

  ## range1 <- ranges[1,] * centerScales.Des[1,2] + centerScales.Des[1,1]
  ## range2 <- ranges[2,] * centerScales.Des[2,2] + centerScales.Des[2,1]

  range1 <- ranges[1,]
  range2 <- ranges[2,]
  
  final.emu.data <- list(mean=emu.mean, var=emu.var, rA=range1, rB=range2, fixed=fixedVals)
  
  implaus <- computeImplaus(obsIndex, final.emu.data, exp.data)

  finalTable <- cbind(pointList, emu.mean, emu.var, implaus)

  
  # we write the table out to a csv file for vis
  write.table(finalTable, file=fname, row.names=FALSE, col.names=c(axisNames, "mean", "variance", "implaus"),
              sep=",", qmethod="double", dec=".")

  invisible(finalTable) # returns final table if you want it
  
}
