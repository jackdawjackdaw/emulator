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
  emu.result <- emulateObsOverDesign(obs, thetas, estim.result$des.scaled,
                                     desA, desB, fixedVal, npts, cov.fn=estim.result$cov.fn,
                                     reg.order=estim.result$reg.order)

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
                        xlabel="", ylabel="", titleIn="", plotDes=FALSE, estim.result, exp.data, feasCut=4){
  npts.side <- 32
  implaus.result <- implausObsOverDesign(obsIndex, dAIndex, dBIndex, fixedVal, estim.result, exp.data, npts=npts.side)
  # cut on some scale (say 5)

  for(i in 1:npts.side){
    for(j in 1:npts.side){
      if(implaus.result$feas[i,j] > feasCut){
        implaus.result$feas[i,j] <- feasCut
      }
    }
  }
    
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

  breaks <- seq(from=0, to=feasCut, length.out=13)
  image(r1, r2, t(implaus.result$feas), axes=FALSE, col=rev(heat.colors(12)), xlab=xlabel, ylab=ylabel)
  contour(r1, r2, t(implaus.result$feas), nlevels=12, col="black", add=TRUE, cex.lab=0.5, labcex=0.8)
  if(plotDes==TRUE){
    points(desA, desB, pch=3)
    title(xlab=xlabel, ylab=ylabel, outer=TRUE, cex.lab=2.0)
    axis(1, cex.axis=1.0)
    axis(2, cex.axis=1.0)
  }
  legend("topright", titleIn, bg="white")
}

##
## like plotImplausOverDesign but we combine all obsIndex and only plot regions where
## the worst  implaus is < feasCut
plotImplausOverDesignCombined <- function(dAIndex, dBIndex, fixedVal, xlabel="", ylabel="",
                                          titleIn="", plotDes=FALSE,
                                          estim.result, exp.data, feasCut=1.8){
  npts.side <- 32
  #browser()
  nobs <- dim(estim.result$train.scaled)[2]
  implaus.list <- vector("list", nobs)

  for(obsIndex in 1:nobs){
    implaus.list[[obsIndex]] <- implausObsOverDesign(obsIndex, dAIndex, dBIndex, fixedVal, estim.result, exp.data, npts=npts.side)
  }

  maxfeas.mat <- matrix(0, nrow=npts.side, ncol=npts.side)
  for(i in 1:npts.side){
    for(j in 1:npts.side){
      maxVal <- 0
      for(obsIndex in 1:nobs){
        if(implaus.list[[obsIndex]]$feas[i,j] > maxVal){
          maxVal <- implaus.list[[obsIndex]]$feas[i,j]
        }
      }
      if(maxVal < feasCut){ 
        maxfeas.mat[i,j] <-  maxVal
      } else {
        maxfeas.mat[i,j] <-  feasCut
      }
    }
  }

  r1 <- implaus.list[[1]]$emu.data$rA
  r2 <- implaus.list[[1]]$emu.data$rB

  desA <- estim.result$des.scaled[,dAIndex]
  desB <- estim.result$des.scaled[,dBIndex]

  desACenterScale <- c(attr(estim.result$des.scaled, "scaled:center")[dAIndex],
                       attr(estim.result$des.scaled, "scaled:scale")[dAIndex])
  
  desBCenterScale <-c(attr(estim.result$des.scaled, "scaled:center")[dBIndex],
                      attr(estim.result$des.scaled, "scaled:scale")[dBIndex])

  desA <- desA * desACenterScale[2] + desACenterScale[1]
  desB <- desB * desBCenterScale[2] + desBCenterScale[1]

  breaks <- seq(from=0, to=feasCut, length.out=13)
  image(r1, r2, t(maxfeas.mat), axes=FALSE, col=rev(heat.colors(12)), xlab=xlabel, ylab=ylabel)
  contour(r1, r2, t(maxfeas.mat), nlevels=8, col="black", add=TRUE, cex.lab=0.5, labcex=0.8)
  if(plotDes==TRUE){
    points(desA, desB, pch=3)
    title(xlab=xlabel, ylab=ylabel, outer=TRUE, cex.lab=2.0)
    axis(1, cex.axis=1.0)
    axis(2, cex.axis=1.0)
  }
  legend("topright", titleIn, bg="white")
  
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

  par(mfrow=c(3,3), mar=c(1,1,0,0), oma=c(4,5,3,0))

  for(i in 0:(nsteps-1)){
    fixV <- minVal + stepSize * i

    # now if we're unscaling we need to unscale the fixed value also
    fixScale <- attr(estim.result$des.scaled, "scaled:scale")[stepDim]
    fixCenter <- attr(estim.result$des.scaled, "scaled:center")[stepDim]
    fixVUnscaled <- fixV * fixScale + fixCenter
    
    buffer <- paste(desNames[stepDim], " : ", round(fixVUnscaled,2), sep="")
    if(i == 6){
      desPlot <- TRUE
      plotImplausOverDesign(obsIndex, plotDimA ,plotDimB, c(fixV, fixedVals),
                      xlabel=desNames[plotDimA], ylabel=desNames[plotDimB], titleIn=buffer,
                      desPlot,  estim.result, exp.data)
    } else {
      desPlot <- FALSE
      plotImplausOverDesign(obsIndex, plotDimA ,plotDimB, c(fixV, fixedVals),
                      xlabel="", ylabel="", titleIn=buffer,
                      desPlot,  estim.result, exp.data)
      
    }
  }
  par(mfrow=c(1,1))
  if(is.null(fixedVals)==FALSE){
    buffer <- paste("obs: ", obsIndex, "fixed: ", fixedVals)
  } else {
    buffer <- paste("obs: ", obsIndex)
  }
  title(main=buffer, outer=TRUE)
}


##
## repeatedly calls plotImplausOverDesignCombined creating a set of 9 plots of obsIndex which span the stepping dimension stepDim
##
## fixedVals -> a vector of additional fixed values (empty unless nparams > 3)
stepPlotDimensionImplausComb <- function(plotDimA, plotDimB, stepDim, fixedVals=NULL, nsteps=9,
                                     estim.result, exp.data){

  minVal <- min(estim.result$des.scaled[,stepDim])
  maxVal <- max(estim.result$des.scaled[,stepDim])
  
  stepSize <- (maxVal - minVal) / nsteps

  par(mfrow=c(3,3), mar=c(1,1,0,0), oma=c(4,5,4,1))

  for(i in 0:(nsteps-1)){
    fixV <- minVal + stepSize * i

    # now if we're unscaling we need to unscale the fixed value also
    fixScale <- attr(estim.result$des.scaled, "scaled:scale")[stepDim]
    fixCenter <- attr(estim.result$des.scaled, "scaled:center")[stepDim]
    fixVUnscaled <- fixV * fixScale + fixCenter
    
    buffer <- paste(desNames[stepDim], " : ", round(fixVUnscaled,2), sep="")
    if(i == 6){

      desPlot <- TRUE
      plotImplausOverDesignCombined(plotDimA ,plotDimB, c(fixV, fixedVals),
                      xlabel=desNames[plotDimA], ylabel=desNames[plotDimB], titleIn=buffer,
                      desPlot,  estim.result, exp.data)
    } else {
      desPlot <- FALSE
      plotImplausOverDesignCombined(plotDimA ,plotDimB, c(fixV, fixedVals),
                      xlabel="", ylabel="", titleIn=buffer,
                      desPlot,  estim.result, exp.data)
      
    }
  }
  par(mfrow=c(1,1))
  title(main="max implausibility < 1.8", outer=TRUE)

}




# I = (EmuMean  - ExpValue)**2 / (Var_Model + Var_Emulator + Var_Data)
# for facundos data we have no var_model
# 
# here we're supposing that expData$obsError is in the form of a standard error
# should support including errors from the model too
computeImplaus <- function(obsIndex, emu.data, expData){
  num <- (emu.data$mean - expData$obsValue[obsIndex])**2

  denom <- (emu.data$var + (expData$obsError[obsIndex])**2)

#  browser()
  
  if(is.null(expData$errModel.sys) == FALSE){
    #we'll model the sys error as normally distributed noise with mean zero
    #the value given in the expData list is taken as a confidence interval (%) around
    #anything we get from the model
    
    sig <- expData$errModel.sys

    # compute the "nvalue" for a confidence interval
    # from the inverse error function (see help(pnorm))
    erfinv <- function(x) { qnorm((1+x)/2)/sqrt(2)}
    nvalue = sqrt(2)*erfinv((1-sig))
    
    sigmaSys = (sig*emu.data$mean) / nvalue
    denom <- denom + (sigmaSys**2)


  }

  # it is not clear how to include the spread or statistical errors in the model output yet
  # this should probably be pushed into the nugget and then the uncertainty will be automatically
  # propagated by the emulator

  # scalar matrix divide i hope?
  I <- sqrt(num / denom)

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
    fixedValSeq <- matrix(0, nrow=nEmuPts, ncol=length(fixedVals))
    for(i in 1:nEmuPts){
      fixedValSeq[i,] <- fixedVals
    }
    pointList <- as.matrix(cbind(pointGrid, fixedValSeq))
  } 
  
  pointList <- as.matrix(pointGrid)


  if(estim.result$cov.fn == 1){
    nthetas <- nParams + 2
  } else {
    nthetas <- 3
  }
 
  
  emuRes <- callEmulateAtList( list(xmodel=des, training=obs),
                              estim.result$thetas[obsIndex,], pointList,
                              nemupts = nEmuPts, nmodelpoints = nModelPoints,
                              nparams = nParams, nthetas = nthetas,
                              cov.fn = estim.result$cov.fn,
                              reg.order = estim.result$reg.order)

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

## compute the implaus over variables A,B,C for each observable in a set of  grids spanning the space
## compute the max of the set of implausibilities for each point in the grid
##
## profit?
gridImplausComb <- function(estim.result=estimResult, exp.data, thresh=1.8, dimA, dimB, dimC, fixedVals=NULL, nGridPts=32,
                            unscale = TRUE, fname="grid-sweep-comb.csv")
{
  nobs <- dim(estim.result$train.scaled)[2]
  ndim <- 3
  imp.list <- vector("list", nobs)
  # this will create temp.csv
  for(obsIndex in 1:nobs){
    # gridImplausSweep produces a matrix which is : pointList, emu.mean, emu.var, implaus
    imp.list[[obsIndex]] <- gridImplausSweep(estim.result=estim.result, exp.data, obsIndex, dimA, dimB, dimC, fixedVals, nGridPts,
                                              unscale=TRUE, fname="temp.csv")
  }

  nparams <- dim(estim.result$des.scaled)[2]
  imp.col <- 3+nparams

  feas.mat <- matrix(0, nrow=nGridPts**ndim, ncol=(nparams+1))
  
  for(i in 1:nGridPts**3){
    maxVal <- 0
    for(obsIndex in 1:nobs){
      if(imp.list[[obsIndex]][i, imp.col] > maxVal){
        maxVal <- imp.list[[obsIndex]][i, imp.col]
      }
    }
    # we do the thresholding in paraview for now
    feas.mat[i,nparams+1] <- maxVal
    feas.mat[i, 1:nparams] <- imp.list[[1]][i, 1:nparams]
  }

  
  dimList <- c(dimA, dimB, dimC)
  axisNames <- rep(NA, ndim)
  dimnames <- attr(estim.result$des.scaled, "dimnames")[[2]]
  for(i in 1:ndim){
    axisNames[i] <- dimnames[dimList[i]]
  }
  
  # we write the table out to a csv file for vis
  write.table(feas.mat, file=fname, row.names=FALSE, col.names=c(axisNames, "max-implaus"),
              sep=",", qmethod="double", dec=".")

}
