source("sa-fns.R")

doSA <- function(nDimensions, nModelPts, nRegresionFns, designMatrix, trainingVector, thetas){
  inputSampleLength <- 100
                                        # the parameters are taken to be all N(0,1) 
  paramVariances <- diag(nDimensions)
  gpEmuVariances <- diag(nDimensions)
  diag(gpEmuVariances) <- thetas[3:(nDimensions+2)]

  x <- seq(-2.5, 2.5, length=inputSampleLength)
  posteriorMean <- rep(0, 100)
  effect <- rep(0, 100)

  # emulator fns
  cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
  cinv <- solve(cMatrix)
  hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
  betaVec <- makeBetaVector(hmatrix, cinv, training)

  r <- rZero(nreg)
  t <- tZero(designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
  emat <- makeEVector(betaVec, hmatrix, training, cinv)

  resultsMeans <- matrix(0, ncol=inputSampleLength, nrow=nDimensions)
  resultsEffects <- matrix(0, ncol=inputSampleLength, nrow=nDimensions)

  for(j in 1:nDimensions){
    for(i in 1:inputSampleLength){
      r1 <- rpOne(x[i], j, nreg)
      t1 <- tpOne(x[i], j, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
      posteriorMean[i] <- t1 %*% emat + r1 %*% betaVec
      effect[i] <- (r1 - r)%*% betaVec + (t1-t) %*% emat
    }
    resultsMeans[j,] <- posteriorMean
    resultsEffects[j,] <- effect
  }
  #resultsMeans
  resultsEffects
}

  


## setup the constants for our particular run
nDimensions <- 6
nModelPts <- 250
nreg <- nDimensions + 1
nThetas <- nDimensions + 2


# now create arrays etc
thetas <- rep(0, nThetas)
designMatrix <- as.matrix(read.table('6d-oak-ohagan-data-sample.txt')[,1:(nDimensions)])
training <- read.table('6d-oak-ohagan-data-sample.txt')[,(nDimensions+1)]
# the parameters are taken to be all N(0,1) 
paramVariances <- diag(nDimensions)
# load the thetas
thetaTemp <- t(read.table('thetas-6.txt'))
for(i in 1:nThetas){
  if(i != 2){
    thetas[i] <- exp(thetaTemp[i])
  } else {
    # nugget is not rescaled, annooooying
    thetas[i] <- thetaTemp[i]
  }
}
x <- seq(-2.5, 2.5, length=100)


plotEffects <- function(x, resultsEffects, ymax=10, ymin=-10){
  Plot(x, resultsEffects[1,], type="l", lty=1, ylim=c(ymin,ymax))
  for(i in 2:5)
    lines(x, resultsEffects[i,], lty=1)
  for(i in 6:10)
    lines(x, resultsEffects[i,], lty=2, col="blue")
  for(i in 11:15)
    lines(x, resultsEffects[i,], lty=3, col="red")
}

plotMeans <- function(x, resultsMeans, ymax=10, ymin=-10){
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
  lines(x, resultsEffects[4,], lty=3, col="red")
  lines(x, resultsEffects[5,], lty=2, col="blue")
  lines(x, resultsEffects[6,], lty=3, col="red")
  
}

