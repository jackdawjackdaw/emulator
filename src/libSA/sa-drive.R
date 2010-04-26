source("sa-fns.R")

nDimensions <- 15
nModelPts <- 250
nreg <- nDimensions + 1

x <- seq(-2.5, 2.5, length=100)

designMatrix <- as.matrix(read.table('15d-oak-ohagan-data-sample.txt')[,1:15])
training <- read.table('15d-oak-ohagan-data-sample.txt')[,16]
# the parameters are taken to be all N(0,1) 
paramVariances <- diag(ndesign)
# load the thetas
thetas <- read.table('thetas-15d.txt')$V1

gpEmuVariances <- diag(ndesign)
diag(gpEmuVariances) <- thetas

# support fns
cMatrix <- makeCovMatrix(designMatrix, thetas, nModelPts)
cinv <- solve(cMatrix)
hmatrix <- makeHMatrix(designMatrix, nModelPts, nreg)
betaVec <- makeBetaVector(hmatrix, cinv, training)

# do p =1 to test 
r <- rZero(nreg)
r1 <- rpOne(x, 1, nreg)
t1 <- tpOne(x, 1, designMatrix, nModelPts, paramVariances, gpEmuVariances)
emat <- makeEVector(beta, hmatrix, training, cinv)


posteriorMean <- r1 %*% beta  + t1 %*% emat

