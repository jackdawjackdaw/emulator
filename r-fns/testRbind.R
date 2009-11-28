source("EmuRbind.R")

## some test functions and the interface to rbind (at the end of the file)
## the main ones of interest are callcode and setEmulatorOptions
##
## setEmulatorOptions -> lets you change the covariance function and lets you change the alpha in the cov fn
## you have to call this sometime before you run callcode! if you don't want to change things you should use
## setDefaultOptions (as in the examples)
##
## callcode -> takes a data frame for the model (and some other stuff) and produces another data frame with the emulated
## mean and variance

# this is a simple example of how to use the callcode function
# we load some data from our general test file and then call the
# rbind lib with the default params
testIsing <- function(){
  model <- readIsingData()
  nmodelpts <- dim(model)[1]
  nemupts <- 100
  setDefaultOps()
  ans<-callcode(model, nmodelpts, nemupts=nemupts)
  plotResults(model, ans)  
  #ans
}


## plots an array of grophs, toprow using lhs sampled input
## bottom row using uniformly sampled input
##
## interstingly (uninterestingly) the uniformly sampled input
## really really sucks
makeSuperPlot <- function(){
  par(mfrow=c(2,5))
  #postscript("superPlot.ps")
  for(i in 1:5)
    testNModelPts(i+3,1)
  for(i in 1:5)
    testNModelPts(i+3,0)
  #dev.off()
}



## change the alpha value for the gauss covFn
## keep the number of points fixed,
## should use a single lhs sample for all the pulls but
varyAlpha <- function(){
  par(mfrow=c(1,5))
  title <- "Gauss Plot"
  for(i in 1:5){
    alpha <- i/5 + 1
    setEmulatorOptions(1, alpha)
    string <- toString(alpha)
    testNModelPtsGauss(8, lhs=0, string)
  }
}


## do a row wise compare of
## gaussian and matern cov fns for different m
compareCovFns <- function(){
  par(mfrow=c(2,5))
  for(i in 1:5)
  # use the gaussian cov fn
    testNModelPtsGauss(i+3, lhs=1, "Gaussian Covariance Fn")

  for(i in 1:5)
  # now set the matern cov fn
    testNmodelPtsMatern(i+3, lhs=1, "Matern Covariance Fn")
}

# test the linear interpolation model,
# it works ok but this is not a good demo yet
testInt <- function(m){
  model <- demoModel(m, 1)
  f <- lm(model$training ~ model$xmodel + I(model$xmodel^2) + I(model$xmodel^3) )
  plot(model, ylim=range(0,6))
  #print(summary(f) )
  #abline(f)
  x <- seq(0,1,by=0.01)
  y<-rep(0, 101)
  y<-f[1]$coefficients[1]*x^0 + f[1]$coefficients[2]*x + f[1]$coefficients[3]*(x^2) + f[1]$coefficients[4]*x^3
  #browser()
  lines(x,y)
  ModelPlus <- demoModel(m+2, 1)
  points(ModelPlus, ylim=range(0,6), col="red")
}


## interpolate (doesnt' work)
## # y is a vec of values of the fn (evaluated at the pts x)
## interpolate <- function(y,x, order, npts, rangemin=0.0, rangemax=1.0){
##   d <- length(x)
##   F = matrix(0, d, order)
##   for(i in 1:d){
##     for(j in order:0){
##       #print(j)
##       F[i,order-j] <- x[i]^(j)
##     }
##   }
##   ##this gives you ans where: F %*% ans = y 
##   ans <- solve(F, y)
##   ##print(y)
##   ##print(F %*% ans)
##   ## fullX <- rep(0, npts)
##   ## for(i in 1:npts)
##   ##   fullX[i] <- i*(1/npts)
##   ## fullY <- rep(0, npts)
##   ## temp <- rep(0, order)
##   ## #browser()
##   ## ## now actually use the thing
##   ## for(i in 1:npts){
##   ##   fullY[i] <- intFn(fullX[i], ans)
##   ## }
##   ## #browser()
##   ## final <- data.frame(x=fullX, y=fullY)
## }

## intFn <- function(x,a){
##   order <- length(a)
##   temp <- rep(0, order)
##   for(i in 1:order)
##     temp[i] <- (x^(i-1))*a[i]
##   print(temp)n
##   ans <- sum(temp)
## }



## for making some demo plots
yM <- function(z) {5*exp(-3*z)*(sin(z*10)) + 2}

demoModel <- function(m, lhs=0){
  if(lhs == 1){
    params <- maximinLHS(m, 1)
  } else {
    params <- matrix(0,m,1)
    params[1,] <- 0
    for(i in 1:(m-1))
      params[i+1,] <- (1/m)*i
  }
  ymodel = rep(NA, m)
  for(i in 1:m)
    ymodel[i] <- yM(params[i,])
  model <- data.frame(xmodel = params, training = ymodel)
  model
}

testNModelPtsGauss <- function(m, lhs=0, title="test"){
  model <- demoModel(m, lhs=lhs)
  nmodelpts <- m
  nemupts <- 200
  print(model)
  setEmulatorOptions(1,1.9)
  ans <- callcode(model, nmodelpts, nemupts=nemupts, rangemin=0.0, rangemax=1.0)
  sequence <- seq(0.0,1.0, length=nemupts)
  actual <- data.frame(x=sequence, y=yM(sequence))
  plotResultsTest(model, ans, actual, title)

}


testNModelPtsMatern <- function(m, lhs=0, title="test"){
  model <- demoModel(m, lhs=lhs)
  nmodelpts <- m
  nemupts <- 200
  print(model)
  setEmulatorOptions(0,1.9)  
  ans <- callcode(model, nmodelpts, nemupts=nemupts, rangemin=0.0, rangemax=1.0)
  sequence <- seq(0.0,1.0, length=nemupts)
  actual <- data.frame(x=sequence, y=yM(sequence))
  plotResultsTest(model, ans, actual, title)
}




# super super slow(if you use matern)
# same as plot results but now we have an analytic actual model to plot
plotResultsTest <- function(model, results, actual, title="testing"){
  # this works in 1d now
  plot(model$xmodel, model$training, ylim=range(0.0,6.0), xlab="x", ylab="y", pch=19)
  title(main=title)
  lines(actual$x, actual$y, col="green", lwd=2, lty=2)
  lines(results$emulatedx, results$emulatedy, col="red", lwd=2)
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- 0.5*sqrt(results$emulatedvar[i])*1.65585
 
  lines(results$emulatedx, results$emulatedy + confidence, col="red", lty=2)
  lines(results$emulatedx, results$emulatedy - confidence, col="red", lty=2)
  grid()
}

      
plotResults <- function(model, results){
  # this works in 1d now
  plot(model$xmodel, model$training)
  lines(results$emulatedx, results$emulatedy, col="red")
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- 0.5*sqrt(results$emulatedvar[i])*1.65585
 
  lines(results$emulatedx, results$emulatedy + confidence, col="blue")
  lines(results$emulatedx, results$emulatedy - confidence, col="blue")
  grid()
}


yGauss <- function(x,y,x0,y0, a,b,c){
  f <- a*(x-x0)^2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)^2
  return(exp(-f))
}


## these gauss functions don't work yet
## a nice plot can be made : scatterplot3d(demoGauss(20), highlight.3d=TRUE, pch=20)
## but you can't rotate the little guy :(
demoGauss <- function(m){
  # this is a 2d func to we'll use a LHS in 2d
  params <- maximinLHS(m,2)
  result <- rep(0,m)
  for(i in 1:m){
    result[i] <- yGauss(params[i,1], params[i,2],0.5,0.5,100,0.0,100) +
      yGauss(params[i,1], params[i,2],0.1,0.1,75,0,75)
  }
  model <- data.frame(xmodel=params, training=result)
  model
}

## this works now but it doesn't seem to produce the correct output
## rbind is probably not set up for 2d output yet
runGauss <-function(m){
  nmodelpts <- m
  model <- demoGauss(nmodelpts)
  nemupts <- 200
  nthetas <- 5
  setDefaultOps()
  ans <- callcode(model, nmodelpts, nemupts=200, nparams=2, nthetas=nthetas)
  ans
}




test1dSin <- function(){
  nmodelpts <- 40
  model <- test1dRCall(nmodelpts)
  setDefaultOps()
  ans<-callcode(model, nmodelpts, rangemin=0.0, rangemax=60, nemupts=100)
  plotResults(model, ans)
}

readIsingData <- function(){
  temp <- read.table("../model-cut.dat")
  model <- data.frame(xmodel=temp[,1], training=temp[,2])
  model
}

testLikely <- function(){
  vertex <- c(0.3, 0.2, 0.2, 0.0001)
  model <- readIsingData()
  nmodelpts <- dim(model)[1]
  l<-callEvalLikelyhood(model, nmodelpts, vertex)
  l
}

