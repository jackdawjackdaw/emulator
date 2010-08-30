source("EmuRbind.R")

## some test functions and the interface to rbind (at the end of the file)
## the main ones of interest are callcode and setEmulatorOptions
##
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
  ans<-callcode(model, nmodelpts, nemupts=nemupts)
  plotResults(model, ans)  
  #ans
}



## for making some demo plots
yM <- function(z) {5*exp(-3*z)*(sin(z*10)) + 2}

## the data set for this
demoModel <- function(m, lhs=0, rangeMin=0.0, rangeMax=1.0, modelFunc=yM){
  if(lhs == 1){
    params <- (maximinLHS(m, 1))*(rangeMax-rangeMin) + rangeMin
  } else {
    params <- matrix(0,m,1)
    x<-  seq(from=rangeMin,to=rangeMax,length.out=m)
    for(i in 1:m)
      params[i,] <- x[i]

  }
  ymodel = rep(NA, m)
  for(i in 1:m)
    ymodel[i] <- modelFunc(params[i,])
  model <- data.frame(xmodel = params, training = ymodel)
  model
}

testNModelPtsGauss <- function(m, lhs=0, title="test", noSet=0, rangeMin=0.0, rangeMax=1.5, printL=0){
  model <- demoModel(m, lhs=lhs, rangeMin, rangeMax)
  nmodelpts <- m
  nemupts <- 200
  print(model)

  if(printL == 1){
    thetas <- callEstimate(model, nmodelpts)
    ans <- callEmulate(model, thetas, nmodelpts, nemupts=nemupts, rangemin=rangeMin, rangemax=rangeMax)
    
    } else { 
      ans <- callcode(model, nmodelpts, nemupts=nemupts, rangemin=rangeMin, rangemax=rangeMax)
    }
  sequence <- seq(0.0,1.0, length=nemupts)
  actual <- data.frame(x=sequence, y=yM(sequence))
  plotResultsTest(model, ans, actual, title)

  if(printL==1){
    mtext(toString(thetas[4], width=9))
  }
    
}

# same as plot results but now we have an analytic actual model to plot
plotResultsTest <- function(model, results, actual, title="testing"){
  # this works in 1d now
  plot(model$xmodel, model$training, ylim=range(0.0,6.0), xlab="x", ylab="y", pch=19, xlim=range(0.0,1.5))
  title(main=title)
  lines(actual$x, actual$y, col="darkolivegreen", lwd=2, lty=2)
  lines(results$emulatedx, results$emulatedy, col="red", lwd=2)
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- sqrt(results$emulatedvar[i])*1.65585
 
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
  temp <- read.table("~/Projects/Emulator/model-cut.dat")
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

