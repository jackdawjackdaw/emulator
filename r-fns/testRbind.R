dyn.load("../src/libRbind.so")
source("emulator-test-data.R")

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

plotResults <- function(model, results){
  # this works in 1d now
  plot(model$xmodel, model$training)
  lines(results$emulatedx, results$emulatedy, col="red")
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- 0.5*sqrt(results$emulatedvar[i])*1.65585
 
  lines(results$emulatedx, results$emulatedy + confidence, col="blue")
  lines(results$emulatedx, results$emulatedy - confidence, col="blue")  
}

test1dSin <- function(){
  nmodelpts <- 40
  model <- test1dRCall(nmodelpts)
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


# i don't quite understand how to push the results together into a data
# frame i should check on this.
# note that the rangemin/max create a square domain in 2d.
# not sure how to grab the results?
callcode <- function(model, nmodelpts, nparams=1, nthetas=4, nemupts=50, rangemin=0.0, rangemax=4.0){

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
