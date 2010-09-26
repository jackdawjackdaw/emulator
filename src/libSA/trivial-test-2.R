source("sa-fns.R")
nDimensions <- 4
nModelPts <- 60
nSamples <- 64 # how many places to eval the sensitivity fns at
nThetas <- nDimensions + 2
nreg <- 1

thetas <- rep(0, nThetas)
temp <- read.table('trivial-model-2.txt')
designMatrix <- as.matrix(temp[,1:nDimensions])
training <- as.matrix(temp[,(nDimensions+1)])

thetas <- exp(t(read.table('thetas-trivial-2.dat')))

x <- seq(-2.5,2.5, length=nSamples)

paramVariances <- diag(nDimensions)
gpEmuVariances <- diag(nDimensions)
                                        # omg! needs to be like this
diag(gpEmuVariances) <- 1/(thetas[3:(nDimensions+2)]^2)
x <- seq(-2.5, 2.5, length=nSamples)
posteriorMean <- rep(0, nSamples)
effect <- rep(0, nSamples)
variance <- rep(0, nSamples)

                                        # emulator fns
cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
cinv <- solve(cMatrix)
hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
betaVec <- makeBetaVector(hmatrix, cinv, training)
wmatrix <- solve(t(hmatrix) %*% cinv %*% hmatrix)


r <- rZero(nreg)
                                        # need to include the scale too
t <- tZero(designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
emat <- makeEVector(betaVec, hmatrix, training, cinv)

resultsMeans <- matrix(0, ncol=nSamples, nrow=nDimensions)
resultsVariances <- matrix(0, ncol=nSamples, nrow=nDimensions)
resultsEffects <- matrix(0, ncol=nSamples, nrow=nDimensions)

for(j in 1:nDimensions){
  r1 <- rep(0, nreg)
  t1 <- rep(0, nModelPts)
  for(i in 1:nSamples){
    r1 <- rpOne(x[i], j, nreg)
    t1 <- thetas[1]*tpOne(x[i], j, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
    posteriorMean[i] <- t1 %*% emat  +  r1 %*% betaVec
    effect[i] <- (r1 - r)%*% betaVec + (t1-t) %*% emat
    variance[i] <- (thetas[1]*upq(x[i], x[i], j, j, nDimensions, paramVariances, gpEmuVariances) - t1 %*% cinv %*% t1 + (r1 - t1 %*% cinv %*% hmatrix) %*% wmatrix %*% t(r1 - t1 %*% cinv %*% hmatrix))
  }
  #browser()
  ## the spinning top, made a sound, like a train across the valley
  resultsVariances[j,] <- variance
  resultsMeans[j,] <- posteriorMean
  resultsEffects[j,] <- effect
  posteriorMean <- rep(0, nSamples)
  variance <- rep(0, nSamples)
}


## do the analytic calcg
an1 <- function(x){
  exp(-0.1*x)*exp(-12.5) + 0.2
}

an2 <- function(x){
  cos(5*x)*exp((-0.1)^2/2) + 0.2
}

an3 <- function(x){
  exp(-12.5)*exp((-0.1)^2/2) + 3*sin(x) + 3*x + 0.2
}

an4 <- function(x){
  exp(-12.5)*exp((-0.1)^2/2) + 0.2*x^2
}

pdf("sensitivity-test.pdf")
plot(x, resultsMeans[1,], type="l", col="black", ylim=c(-5,5), xlab="standardised input value", ylab="E(Y|x_i)",
     main = "Sensitivity Analysis of y = exp(-0.1x)cos(5y) + 3sin(z) + 3z + 0.2w^2")
lines(x, resultsMeans[2,], col="red")
lines(x, resultsMeans[3,], col="green")
lines(x, resultsMeans[4,], col="purple")

points(x, an1(x), col="black")
points(x, an2(x), col="red")
points(x, an3(x), col="green")
points(x, an4(x), col="purple")

abconf <- sqrt(abs(resultsVariances))*0.5*1.65585
lines(x, resultsMeans[1,] + abconf[1,], col="black", lty=2)
lines(x, resultsMeans[1,] - abconf[1,], col="black", lty=2)
lines(x, resultsMeans[2,] + abconf[2,], col="red", lty=2)
lines(x, resultsMeans[2,] - abconf[2,], col="red", lty=2)
lines(x, resultsMeans[3,] + abconf[3,], col="green", lty=2)
lines(x, resultsMeans[3,] - abconf[3,], col="green", lty=2)
lines(x, resultsMeans[4,] + abconf[4,], col="purple", lty=2)
lines(x, resultsMeans[4,] - abconf[4,], col="purple", lty=2)


legend(x=1,y=-1, legend=c('E*{E(Y|x)}', 'E*{E(Y|y)}', 'E*{E(Y|z)}', 'E*{E(Y|w)}', 'E(Y|x)', 'E(Y|y)', 'E(Y|z)', 'E(Y|w)'),
       col=c('black', 'red', 'green', 'purple', 'black', 'red', 'green', 'purple'),
       pch=c(-1,-1,-1,-1,1,1,1,1),
       lty=c(1,1,1,1,0,0,0,0),
       bg='white'
       )
dev.off()       
       
