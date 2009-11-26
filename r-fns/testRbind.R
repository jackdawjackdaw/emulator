dyn.load("../src/libRbind.so")
source("emulator-test-data.R")
library("lhs")
library("scatterplot3d")

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
  #postscript("psuperPlot.ps")
  for(i in 1:5)
    testNModelPts(i+3,1)
  for(i in 1:5)
    testNModelPts(i+3,0)
  #dev.off()
}


compareCovFns <- function(){
  par(mfrow=c(1,2))
  m <- 8
  # use the gaussian cov fn
  testNModelPts(m, lhs=1)
  # now set the matern cov fn
  setEmulatorOptions(0,1.9)
  model <- demoModel(m, lhs=1)
  nmodelpts <- m
  nemupts <-  200
  ansMatern <- callcode(model, nmodelpts, nemupts=nemupts, rangemin=0.0, rangemax=1.0)
  sequence <- seq(0.0,1.0, length=nemupts)
  actual <- data.frame(x=sequence, y=yM(sequence))
  plotResultsTest(model, ansMatern, actual)
}
                       



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

testNModelPts <- function(m, lhs=0){
  model <- demoModel(m, lhs=lhs)
  nmodelpts <- m
  nemupts <- 200
  print(model)
  setDefaultOps()
  ans <- callcode(model, nmodelpts, nemupts=nemupts, rangemin=0.0, rangemax=1.0)
  sequence <- seq(0.0,1.0, length=nemupts)
  actual <- data.frame(x=sequence, y=yM(sequence))
  plotResultsTest(model, ans, actual)

}




# super super slow(if you use matern)
# same as plot results but now we have an analytic actual model to plot
plotResultsTest <- function(model, results, actual){
  # this works in 1d now
  plot(model$xmodel, model$training, ylim=range(0.0,6.0))
  lines(actual$x, actual$y, col="green")
  lines(results$emulatedx, results$emulatedy, col="red")
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- 0.5*sqrt(results$emulatedvar[i])*1.65585
 
  lines(results$emulatedx, results$emulatedy + confidence, col="blue")
  lines(results$emulatedx, results$emulatedy - confidence, col="blue")
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

runGauss <-function(m){
  nmodelpts <- m
  model <- demoGauss(nmodelpts)
  nemupts <- 200
  nthetas <- 5
  ans <- callcode(model, nmodelpts, nemupts=200, nparams=2, nthetas=nthetas)
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

setDefaultOps <- function(){
  setEmulatorOptions(1,1.9)
}


## i've taken away the call in the C function, so you need to do this
## before you CALL CODE or DEATH
# the function we're calling has no return, it just switches
# some options around so this should be pretty simple
setEmulatorOptions <- function(useGauss, alpha){
  foo<-.C("setEmulatorOptions",
     as.integer(useGauss),
     as.double(alpha))

}


# i don't quite understand how to push the results together into a data
# frame i should check on this.
# note that the rangemin/max create a square domain in 2d.
# not sure how to grab the results?
callcode <- function(model, nmodelpts, nparams=1, nthetas=4, nemupts=50, rangemin=0.0, rangemax=4.0){

  if(nparams==1){
  
    res<-  .C("callEmulator",
     as.double(model$xmodel),
     as.integer(nparams),
     as.double(model$training),
     as.integer(nmodelpts),
     as.integer(nthetas),
     finalx = double(nemupts*nemupts),
     as.integer(nemupts),
     finaly = double(nemupts),
     finalvar = double(nemupts),
     as.double(rangemin),
     as.double(rangemax))
  }else if(nparams==2){
    newmodel <- rep(NA, 2*nmodelpts)
    print(model)
    # interleave the xmodel array
    for(i in 1:nmodelpts)
      newmodel[2*i-1] <- model$xmodel.1[i]
    for(i in 1:nmodelpts)
      newmodel[2*i] <- model$xmodel.2[i]
    print(newmodel)
    
    res<-  .C("callEmulator",
     as.double(newmodel),
     as.integer(nparams),
     as.double(model$training),
     as.integer(nmodelpts),
     as.integer(nthetas),
     finalx = double(nemupts*nemupts),
     as.integer(nemupts),
     finaly = double(nemupts),
     finalvar = double(nemupts),
     as.double(rangemin),
     as.double(rangemax))
    

  } else {
    print("sorry, won't work with nparams > 2")
  }

  results <- data.frame(emulatedx=res$finalx[1:nemupts], emulatedy=res$finaly, emulatedvar=res$finalvar)
  results
} 

## just estimates the thetas for a model (this is the slow ass part)
callEstimate <- function(model, nmodelpts,nparams=1, nthetas=4){

  res <- .C("callEstimate",
            as.double(model$xmodel),
            as.integer(nparams),
            as.double(model$training),
            as.integer(nmodelpts),
            as.integer(nthetas),
            thetas = double(nthetas))
  res$thetas
}

## test this, its not working right
testCallEm <- function(){
  m <- 6
  model <- demoModel(6, 0)
  # just made up but about right for matern
  ans <- c(0.89, 0.54, 0.334, 0.64)
  f1<-callEmulate(model, ans, m, nemupts=10)
  print(f1)
  f2<-callEmulate(model, ans,m, nemupts=10)
  print(f2)
}

## use a given set of thetas to emulate the code
callEmulate <- function(model, thetas, nmodelpts, nparams=1, nthetas=4, nemupts=100, rangemin=0.0, rangemax=1.0){
  res <- .C("callEmulate",
            as.double(model$xmodel),
            as.integer(nparams),
            as.double(model$training),
            as.integer(nmodelpts),
            as.double(thetas),
            as.integer(nthetas),
            finalx = double(nemupts*nemupts),
            as.integer(nemupts),
            finaly = double(nemupts),
            finalvar = double(nemupts),
            as.double(rangemin),
            as.double(rangemax))
            
   results <- data.frame(emulatedx=res$finalx[1:nemupts], emulatedy=res$finaly, emulatedvar=res$finalvar)           
  results
}

   
callEvalLikelyhood <- function(model, nmodelpoints, vertex, nparams=1,nthetas=4){
  answer <- 0.0
  likely <- .C("callEvalLikelyhood",
               as.double(model$xmodel),
               as.integer(nparams),
               as.double(model$training),
               as.integer(nmodelpoints),
               as.integer(nthetas),
               as.double(vertex),
               as.double(answer))
  likely[[7]]
}
