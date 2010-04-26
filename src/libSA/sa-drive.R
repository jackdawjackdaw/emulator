source("sa-fns.R")

nDimensions <- 15
nModelPts <- 250
nreg <- nDimensions + 1

x <- seq(-2.5, 2.5, length=100)
posteriorMean <- rep(0, 100)

designMatrix <- as.matrix(read.table('15d-oak-ohagan-data-sample.txt')[,1:15])
training <- read.table('15d-oak-ohagan-data-sample.txt')[,16]
# the parameters are taken to be all N(0,1) 
paramVariances <- diag(nDimensions)
# load the thetas
thetas <- t(read.table('thetas.txt')[3:17])
# cheat
#thetas <- rep(1.0, 15)

gpEmuVariances <- diag(nDimensions)
diag(gpEmuVariances) <- thetas

# support fns
cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
cinv <- solve(cMatrix)
hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
betaVec <- makeBetaVector(hmatrix, cinv, training)

# do p =1 to test 
r <- rZero(nreg)
r1 <- rpOne(0.2, 1, nreg)
t1 <- tpOne(0.3, 1, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
emat <- makeEVector(betaVec, hmatrix, training, cinv)

for(i in 1:100){
  r1 <- rpOne(x[i], 2, nreg)
  t1 <- tpOne(x[i], 2, designMatrix, nModelPts, nDimensions, paramVariances, gpEmuVariances)
  posteriorMean[i] <- t1 %*% emat #+ r1 %*% betaVec
}



