source("sa-fns.R")
nDimensions <- 2
nModelPts <- 30
nSamples <- 64 # how many places to eval the sensitivity fns at
nThetas <- nDimensions + 2
nreg <- 1

thetas <- rep(0, nThetas)
temp <- read.table('trivial-model.txt')
designMatrix <- as.matrix(temp[,1:nDimensions])
training <- as.matrix(temp[,(nDimensions+1)])

thetas <- exp(t(read.table('thetas-trivial.txt')))

x <- seq(-2.5,2.5, length=nSamples)

paramVariances <- diag(nDimensions)
gpEmuVariances <- diag(nDimensions)
                                        # omg! needs to be like this
diag(gpEmuVariances) <- 1/(thetas[3:(nDimensions+2)]^2)
x <- seq(-2.5, 2.5, length=nSamples)
posteriorMean <- rep(0, nSamples)
effect <- rep(0, nSamples)

                                        # emulator fns
cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
cinv <- solve(cMatrix)
hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
betaVec <- makeBetaVector(hmatrix, cinv, training)

r <- rZero(nreg)
                                        # need to include the scale too
t <- tZero(designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
emat <- makeEVector(betaVec, hmatrix, training, cinv)

resultsMeans <- matrix(0, ncol=nSamples, nrow=nDimensions)
resultsEffects <- matrix(0, ncol=nSamples, nrow=nDimensions)

for(j in 1:nDimensions){
  r1 <- rep(0, nreg)
  t1 <- rep(0, nModelPts)
  for(i in 1:nSamples){
    r1 <- rpOne(x[i], j, nreg)
    t1 <- thetas[1]*tpOne(x[i], j, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
    posteriorMean[i] <- t1 %*% emat  +  r1 %*% betaVec
    effect[i] <- (r1 - r)%*% betaVec + (t1-t) %*% emat
  }
  resultsMeans[j,] <- posteriorMean
  resultsEffects[j,] <- effect
  posteriorMean <- rep(0, 64)
}


## do the analytic calc
an1 <- function(x){
  0.1*cos(x)+0.4
}

an2 <- function(x){
  0.4*x^2 + 0.1*exp(-1/2)
}


plot(x, resultsMeans[1,], type="l", col="black", ylim=c(0,2))
lines(x, resultsMeans[2,], col="red")
points(x, an1(x), col="black")
points(x, an2(x), col="red")


