
# the data frame selection forces these to come in as doubles
a1 <- read.table('a1-vector.dat')$V1
a2 <- read.table('a2-vector.dat')$V1
a3 <- read.table('a3-vector.dat')$V1
m  <- as.matrix(read.table('m-matrix.dat'))

# set the dimension
ndim <- 15

x <- seq(-2.5, 2.5, length=256)

evalue <- function(x,p, a1,a2,a3,m){
  part1 <- a1[p]*x
  part2 <- a2[p]*sin(x)
  part3 <- a3[p]*cos(x)
  temp <- rep(0, length(a3)-1)
  for(i in 1:length(a3)){
    if(i != p){
      temp[i] <- a3[i]*exp(-1/2)
    }
  }
  part31 <- sum(temp)
  part4 <- x*m[p,p]*x
  for(i in 1:length(a3)){
    temp[i] <- 0.0
    if(i != p){
      temp[i] <- m[i,i]
    }
  }
  part41 <- sum(temp)
  part1+part2+part3+part31+part4+part41
}        
  
resultsAna <- matrix(0, ncol=256, nrow=15)
for(j in 1:15){
    resultsAna[j,] <- evalue(x,j,a1,a2,a3,m)
}

plot(x,resultsAna[15,], type="l", col="red", xlab="x", ylab="E(y|x)")
for(j in 14:1){
  if(j > 10){
    lines(x, resultsAna[j,], col="red")
  } else if(j > 5) {
    lines(x, resultsAna[j,], col="blue")
  } else {
    lines(x, resultsAna[j,], col="black")
  }
}
