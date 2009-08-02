dyn.load("../libRbind.so")

# this is a simple example of how to use the callcode function
# we load some data from our general test file and then call the
# rbind lib with the default params
testIsing <- function(){
  model <- read.table("../../model-cut.dat")
  nmodelpts <- dim(model)[1]
  nemupts <- 100
  ans<-callcode(model, nmodelpts, nemupts=nemupts)
  lines(ans$finalx[1:nemupts], ans$finaly, col="red")
}

test1dSin <- function(){
  nmodelpts <- 50
  model <- test1dRCall(nmodelpts)
  plot(model)
  ans<-callcode(model, nmodelpts, rangemin=0.0, rangemax=60)
  # there seems to be a silly bug in the code here
  lines(ans$finalx[1:50], ans$finaly, col="red")
}



# i don't quite understand how to push the results together into a data
# frame i should check on this.
# note that the rangemin/max create a square domain in 2d.
callcode <- function(model, nmodelpts, nparams=1, nthetas=4, nemupts=50, rangemin=0.0, rangemax=4.0){

  .C("callEmulator",
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

  #results <- data.frame(emulatedx=finalx, emulatedy=finaly, emulatedvar=finalvar)
  #results
} 
   
   
