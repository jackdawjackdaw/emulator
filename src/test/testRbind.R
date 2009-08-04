dyn.load("../libRbind.so")

# this is a simple example of how to use the callcode function
# we load some data from our general test file and then call the
# rbind lib with the default params
testIsing <- function(){
  model <- read.table("../../model-cut.dat")
  nmodelpts <- dim(model)[1]
  nemupts <- 100
  ans<-callcode(model, nmodelpts, nemupts=nemupts)
  plotResults(model, ans)  
  ans
}

plotResults <- function(model, results){
  plot(model)
  lines(results$emulatedx, results$emulatedy, col="red")
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- 0.5*sqrt(results$emulatedvar[i])*1.65585
  lines(results$emulatedx, results$emualtedy + confidence, col="blue")
  lines(results$emulatedx, results$emualtedy - confidence, col="blue")  
}

test1dSin <- function(){
  nmodelpts <- 50
  model <- test1dRCall(nmodelpts)
  ans<-callcode(model, nmodelpts, rangemin=0.0, rangemax=60)
  plotResults(model, ans)
}



# i don't quite understand how to push the results together into a data
# frame i should check on this.
# note that the rangemin/max create a square domain in 2d.
# not sure how to grab the results?
callcode <- function(model, nmodelpts, nparams=1, nthetas=4, nemupts=50, rangemin=0.0, rangemax=4.0){

res<-  .C("callEmulator",
     as.double(model[,1]),
     as.integer(nparams),
     as.double(model[,2]),
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
   
   
