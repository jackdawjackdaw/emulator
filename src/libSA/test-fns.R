source("sa-fns.R")
## testSet <- as.matrix(read.table('mathematica-model.txt'))
## xmodel <- testSet[,1]
## training <- testSet[,2]

## thetaTest <- c(0.007, 0.00383, 0.929)

## cm <- makeCovMatrix(as.matrix(xmodel), thetaTest, length(xmodel))

## hm <- makeHMatrix(as.matrix(xmodel), 19, 2)

## bv <- makeBetaVector(hm, solve(cm), training)

## ev <- makeEVector(bv, hm, training, solve(cm))

## test the t function
nDimensions <- 15
nModelPts <- 250
nreg <- nDimensions + 1
nThetas <- nDimensions + 2
thetas <- rep(0, nThetas)
designMatrix <- as.matrix(read.table('15d-oak-ohagan-data-sample.txt')[,1:(nDimensions)])
training <- read.table('15d-oak-ohagan-data-sample.txt')[,(nDimensions+1)]
                                        # the parameters are taken to be all N(0,1) 
paramVariances <- diag(nDimensions)

cMatrix <- matrix(0, nrow=250, ncol=250)


thetas[1] <- 1
thetas[2] <- 0
thetas[3] <- thetas[4] <- thetas[5] <- thetas[6] <- thetas[7] <- 10
thetas[8] <- thetas[9] <- thetas[10] <- thetas[11] <- thetas[12] <- 1
thetas[13] <- thetas[14] <- thetas[15] <- thetas[16] <- thetas[17] <- 0.1


gpEmuVariances <- diag(nDimensions)
# omg! needs to be like this
diag(gpEmuVariances) <- 1/(thetas[3:(nDimensions+2)]^2)


                                        # emulator fns
cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
cinv <- solve(cMatrix)
hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
betaVec <- makeBetaVector(hmatrix, cinv, training)
emat <- makeEVector(betaVec, hmatrix, training, cinv)

tZero <- tpOne(0, 1, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
rZero <- rpOne(0, 1, nreg)

zeroVal <- rZero %*% betaVec + tZero %*% emat
