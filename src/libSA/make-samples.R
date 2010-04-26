library("lhs")

## produce lhs samples of the Oak & O'hagan data and use it to test the E*{E(Y|xp)} code in effects.c
## don't panic!
## Test function is y=x'Mx+a1'x+a2'sin(x)+a3'cos(x)
## where M is in m-matrix.dat
## and a1, a2, a3 are in a1-vector.dat etc

# the data frame selection forces these to come in as doubles
a1 <- read.table('a1-vector.dat')$V1
a2 <- read.table('a2-vector.dat')$V1
a3 <- read.table('a3-vector.dat')$V1
m  <- as.matrix(read.table('m-matrix.dat'))

# set the dimension
ndim <- 15

# they use 250 points in the paper
sampleSize <- 250 

params  <- maximinLHS(sampleSize, ndim)*5 - 2.5
designVal <- rep(0, sampleSize)

yvalue <- function(x){
  a1 %*% x + a2 %*% (sin(x)) + a3 %*% (cos(x)) + x %*% m %*% x
}


for(i in 1:sampleSize) {
  designVal[i] <- yvalue(params[i,])
}

designFrame <- data.frame(x=params, y=designVal)

write.table(designFrame,"15d-oak-ohagan-data-sample.txt", row.names=FALSE, col.names=FALSE)


