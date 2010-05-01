library("lhs")

## produce lhs samples of the Oak & O'hagan data and use it to test the E*{E(Y|xp)} code in effects.c
## don't panic!
## Test function is y=x'Mx+a1'x+a2'sin(x)+a3'cos(x)
## where M is in m-matrix.dat
## and a1, a2, a3 are in a1-vector.dat etc

## # the data frame selection forces these to come in as doubles
## a1 <- read.table('a1-vector.dat')$V1
## a2 <- read.table('a2-vector.dat')$V1
## a3 <- read.table('a3-vector.dat')$V1
## m  <- as.matrix(read.table('m-matrix.dat'))

## # set the dimension
## ndim <- 15

## # they use 250 points in the paper
## sampleSize <- 250 

## params  <- maximinLHS(sampleSize, ndim)*5 - 2.5
## designVal <- rep(0, sampleSize)

## yvalue <- function(x){
##   a1 %*% x + a2 %*% (sin(x)) + a3 %*% (cos(x)) + x %*% m %*% x
## }


## for(i in 1:sampleSize) {
##   designVal[i] <- yvalue(params[i,])
## }

## designFrame <- data.frame(x=params, y=designVal)

## write.table(designFrame,"15d-oak-ohagan-data-sample.txt", row.names=FALSE, col.names=FALSE)


## ## now make the reduced data set
## ndimred <- 6
## a1red <- c(a1[1], a1[2], a1[6], a1[7], a1[11], a1[12])
## a2red <- c(a2[1], a2[2], a2[6], a2[7], a2[11], a2[12])
## a3red <- c(a3[1], a3[2], a3[6], a3[7], a3[11], a3[12])

## mred <- matrix(0, nrow=6, ncol=6)
## removeVec <- c(-3,-4,-5,-8,-9,-10,-13,-14,-15)
## for(i in 1:ndimred)
##   mred[i,] <- m[i,removeVec]

## yvalred <- function(x){
##   a1red %*% x + a2red %*% (sin(x)) + a3red %*% (cos(x)) + x %*% mred %*% x
## }

## paramsred  <- maximinLHS(sampleSize, ndimred)*5 - 2.5
## yvalsreduced <- rep(0, sampleSize)
## for(i in 1:sampleSize)
##   yvalsreduced[i] <- yvalred(paramsred[i,])


## designFrameRed <- data.frame(x=paramsred, y=yvalsreduced)

## write.table(designFrameRed,"6d-oak-ohagan-data-sample.txt", row.names=FALSE, col.names=FALSE)

## now make a super simple model
ndim <- 4
sampleSize <- 15*ndim

yTrivial <- function(x){
  y <- 0.1*cos(x[1]) + 0.4*x[2]^2
}
yTrivial2 <- function(x){
  y <- exp(-0.1*x[1])*cos(5*x[2]) + 3*sin(x[3]) + 3*x[3] + 0.2*x[4]^2
}

paramsT <- maximinLHS(sampleSize, ndim)*5 - 2.5
ytriv <- rep(0, sampleSize)

for(i in 1:sampleSize)
  ytriv[i] <- yTrivial2(paramsT[i,])

designFrameTriv <- data.frame(x=paramsT, y=ytriv)
write.table(designFrameTriv, "trivial-model-2.txt", row.names=FALSE, col.names=FALSE)

