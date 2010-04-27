source("sa-fns.R")

nDimensions <- 15
nModelPts <- 250
nreg <- nDimensions + 1

x <- seq(-2.5, 2.5, length=100)
posteriorMean <- rep(0, 100)
effect <- rep(0, 100)
thetas <- rep(0, 17)

designMatrix <- as.matrix(read.table('15d-oak-ohagan-data-sample.txt')[,1:15])
training <- read.table('15d-oak-ohagan-data-sample.txt')[,16]
# the parameters are taken to be all N(0,1) 
paramVariances <- diag(nDimensions)
# load the thetas
thetaTemp <- t(read.table('thetas.txt'))
for(i in 1:17){
  if(i != 1){
    thetas[i] <- exp(thetaTemp[i])
  } else {
    # nugget is not rescaled, annooooying
    thetas[i] <- thetaTemp[i]
  }
}

# cheat
#thetas <- rep(1.0, 15)

gpEmuVariances <- diag(nDimensions)
diag(gpEmuVariances) <- thetas[3:17]

# support fns
cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
cinv <- solve(cMatrix)
hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
betaVec <- makeBetaVector(hmatrix, cinv, training)

# do p =1 to test 
r <- rZero(nreg)
t <- tZero(designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
r1 <- rpOne(0.2, 2, nreg)
t1 <- tpOne(0.2, 2, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
emat <- makeEVector(betaVec, hmatrix, training, cinv)

resultsMeans <- matrix(0, ncol=100, nrow=nDimensions)
resultsEffects <- matrix(0, ncol=100, nrow=nDimensions)

for(j in 1:nDimensions){
  for(i in 1:100){
    r1 <- rpOne(x[i], j, nreg)
    t1 <- tpOne(x[i], j, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
    posteriorMean[i] <- t1 %*% emat + r1 %*% betaVec
    effect[i] <- (r1 - r)%*% betaVec + (t1-t) %*% emat
  }
  resultsMeans[j,] <- posteriorMean
  resultsEffects[j,] <- effect
}



plot(x, resultsEffects[1,], type="l", lty=1, ylim=c(-10,10))
for(i in 2:5)
  lines(x, resultsEffects[i,], lty=1)
for(i in 6:10)
  lines(x, resultsEffects[i,], lty=2, col="blue")
for(i in 11:15)
  lines(x, resultsEffects[i,], lty=3, col="red")

plot(x, resultsMeans[1,], type="l", lty=1, ylim=c(-10,10))
for(i in 2:5)
  lines(x, resultsMeans[i,], lty=1)
for(i in 6:10)
  lines(x, resultsMeans[i,], lty=2, col="blue")
for(i in 11:15)
  lines(x, resultsMeans[i,], lty=3, col="red")


  

