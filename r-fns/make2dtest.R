library("lhs")

yGauss <- function(x,y,x0,y0, a,b,c){
  f <- a*(x-x0)^2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)^2
  return(exp(-f))
}


## these gauss functions don't work yet
## a nice plot can be made : scatterplot3d(demoGauss(20), highlight.3d=TRUE, pch=20)
## but you can't rotate the little guy :(
gaussTest <- function(m){
  # this is a 2d func to we'll use a LHS in 2d
  fullSize <- 500
  paramsSmall <- maximinLHS(m,2) # what we'll test the emu with
  paramsFull <- maximinLHS(fullSize,2) # a pretty full sampling
  resultFull <- rep(0,fullSize)
  resultSmall <- rep(0,m)
  for(i in 1:m){
    resultSmall[i] <- 2*yGauss(paramsSmall[i,1], paramsSmall[i,2],0.5,0.5,40,0,40) +
      3*yGauss(paramsSmall[i,1], paramsSmall[i,2],0.4,0.3,50,20,75) + 0.1*rnorm(1)
  }

  for(i in 1:fullSize){
    resultFull[i] <- 2*yGauss(paramsFull[i,1], paramsFull[i,2],0.5,0.5,40,0,40) +
      3*yGauss(paramsFull[i,1], paramsFull[i,2],0.4,0.3,50,20,75) + 0.1*rnorm(1)
  }

  modelSmall <- data.frame(xmodel=paramsSmall, training=resultSmall)
  modelFull <- data.frame(xmodel=paramsFull, training=resultFull)

  write.table(modelSmall, "2d-gauss-emu-sample.txt", row.names=FALSE, col.names=FALSE)
  write.table(modelFull, "2d-gauss-full.txt", row.names=FALSE, col.names=FALSE)
}
