#dyn.load("~/local/lib/libRBIND.so")
#source("emulator-test-data.R")
library("lhs")
#library("scatterplot3d")

## this selects the power exponential function
## with alpha = 1.9
setDefaultOps <- function(){
  
 
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
     finalx = double(nparams*nemupts),
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
     finalx = double(nparams*nemupts),
     as.integer(nemupts),
     finaly = double(nemupts),
     finalvar = double(nemupts),
     as.double(rangemin),
     as.double(rangemax))
    

  } else {
    print("sorry, won't work with nparams > 2")
  }
  #browser()
  results <- data.frame(emulatedx=res$finalx[1:nemupts], emulatedy=res$finaly, emulatedvar=res$finalvar)
  results
} 

## just estimates the thetas for a model (this is the slow ass part)
callEstimate <- function(model, nmodelpts,nparams=1, nthetas=3){

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
callEmulate <- function(model, thetas, nmodelpts, nparams=1, nthetas=3, nemupts=20, rangemin=0.0, rangemax=1.0){
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

## interpolate the data given by xvec and yvec at the point
## xinterp, we need xvec and yvec to be the same length 
##
callInterpolate <- function(xvec, yvec, xinterp){
  if(length(xvec) == length(yvec)){
    yinterp <- .C("lagrange_interp",
                as.double(xvec),
                as.double(yvec),
                as.integer(length(xvec)),
                as.double(xinterp))
                
  } else {
    print("error xvec not same length as yvec")
    yinterp <- 0.0
  }
  yinterp[[4]][1]
}
                                              
