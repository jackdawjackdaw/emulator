source("sa-fns.R")
plotData <- function(designMatrix, training){
  par(mfrow=c(2,3))
  plot(designMatrix[,1], training, type="p")
  for(i in 2:6)
    plot(designMatrix[,i], training, type="p")
}


plotEffects <- function(x, resultsEffects, ymax=10, ymin=-10){
  Plot(x, resultsEffects[1,], type="l", lty=1, ylim=c(ymin,ymax))
  for(i in 2:5)
    lines(x, resultsEffects[i,], lty=1)
  for(i in 6:10)
    lines(x, resultsEffects[i,], lty=2, col="blue")
  for(i in 11:15)
    lines(x, resultsEffects[i,], lty=3, col="red")
}

plotMeans <- function(x, resultsMeans, ymax=10, ymin=-5){
  plot(x, resultsMeans[1,], type="l", lty=1, ylim=c(ymin,ymax))
  for(i in 2:5)
    lines(x, resultsMeans[i,], lty=1)
  for(i in 6:10)
    lines(x, resultsMeans[i,], lty=2, col="blue")
  for(i in 11:15)
    lines(x, resultsMeans[i,], lty=3, col="red")
}

plotMeansRed <- function(x, resultsMeans, ymax=10, ymin=-10){
  plot(x, resultsMeans[1,], type="l", lty=1, ylim=c(ymin,ymax))
  lines(x, resultsMeans[2,], lty=1)
  lines(x, resultsMeans[3,], lty=2, col="blue")
  lines(x, resultsMeans[4,], lty=2, col="blue")
  lines(x, resultsMeans[5,], lty=3, col="red")
  lines(x, resultsMeans[6,], lty=3, col="red")
}

plotEffectsRed <- function(x, resultsEffects, ymax=10, ymin=-10){
  plot(x, resultsEffects[1,], type="l", lty=1, ylim=c(ymin,ymax))
  lines(x, resultsEffects[2,], lty=1)
  lines(x, resultsEffects[3,], lty=2, col="blue")
  lines(x, resultsEffects[4,], lty=2, col="blue")
  lines(x, resultsEffects[5,], lty=3, col="red")
  lines(x, resultsEffects[6,], lty=3, col="red")
  
}



doSA <- function(nDimensions, nModelPts, nRegresionFns, designMatrix, trainingVector, thetas){
  inputSampleLength <- 100
                                        # the parameters are taken to be all N(0,1) 
  paramVariances <- diag(nDimensions)
  gpEmuVariances <- diag(nDimensions)
  # omg! needs to be like this
  diag(gpEmuVariances) <- 1/(thetas[3:(nDimensions+2)]^2)
  x <- seq(-2.5, 2.5, length=inputSampleLength)
  posteriorMean <- rep(0, 100)
  effect <- rep(0, 100)

  # emulator fns
  cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
  cinv <- solve(cMatrix)
  hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
  betaVec <- makeBetaVector(hmatrix, cinv, training)

  r <- rZero(nreg)
  # need to include the scale too
  t <- tZero(designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
  emat <- makeEVector(betaVec, hmatrix, training, cinv)

  resultsMeans <- matrix(0, ncol=inputSampleLength, nrow=nDimensions)
  resultsEffects <- matrix(0, ncol=inputSampleLength, nrow=nDimensions)

  #browser()


   for(j in 1:nDimensions){
     r1 <- rep(0, nreg)
     t1 <- rep(0, inputSampleLength)
     for(i in 1:inputSampleLength){
       r1 <- rpOne(x[i], j, nreg)
       t1 <- tpOne(x[i], j, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
       posteriorMean[i] <- t1 %*% emat  +  r1 %*% betaVec
       effect[i] <- (r1 - r)%*% betaVec + (t1-t) %*% emat
     }
     resultsMeans[j,] <- posteriorMean
     resultsEffects[j,] <- effect
     posteriorMean <- rep(0, 100)
   }
   resultsMeans
   #resultsEffects
 }




 ## setup the constants for our particular run
 nDimensions <- 15
 nModelPts <- 250
 nreg <- nDimensions + 1
 nThetas <- nDimensions + 2


 # now create arrays etc
thetas <- rep(0, nThetas)
temp <- read.table('15d-oak-ohagan-data-sample')
designMatrix <- as.matrix(temp[,1:(nDimensions)])
training <- temp[,(nDimensions+1)]
 # the parameters are taken to be all N(0,1) 
paramVariances <- diag(nDimensions)
 # load the thetas
 ## thetaTemp <- t(read.table('thetas.txt'))
 ## for(i in 1:nThetas){
 ##   if(i != 2){
 ##     thetas[i] <- exp(thetaTemp[i])
 ##   } else {
 ##     # nugget is not rescaled, annooooying
 ##     thetas[i] <- thetaTemp[i]
 ##   }
 ## }

 ## cheat
thetas[1] <- 0.25
thetas[2] <- 0
thetas[3] <- thetas[4] <- thetas[5] <- thetas[6] <- thetas[7] <- 1
thetas[8] <- thetas[9] <- thetas[10] <- thetas[11] <- thetas[12] <- 1
thetas[13] <- thetas[14] <- thetas[15] <- thetas[16] <- thetas[17] <- 1


x <- seq(-2.5, 2.5, length=100)

means <- doSA(nDimensions, nModelPts, nreg, designMatrix, training, thetas)

#postscript("analytic-emulated-2up.ps")
par(mfrow=c(1,2))
plotMeans(x, means,13,0)
source("analytic-sa.R")


