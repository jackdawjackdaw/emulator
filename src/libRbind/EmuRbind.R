library("lhs")

# i don't quite understand how to push the results together into a data
# frame i should check on this.
# note that the rangemin/max create a square domain in 2d.
# not sure how to grab the results?
callcode <- function(model, nmodelpts, nparams=1, nthetas=3, nemupts=50, rangemin=0.0, rangemax=4.0){

  if(nparams==1){
  
    res<-  .C("callEmulator",
     as.double((model$xmodel)),
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
    ## \todo fix emulator to work with n >> 2
    print("sorry, won't work with nparams > 2")
  }
  #browser()
  results <- data.frame(emulatedx=res$finalx[1:nemupts], emulatedy=res$finaly, emulatedvar=res$finalvar)
  results
} 

## just estimates the thetas for a model (this is the slow ass part)
callEstimate <- function(model, nmodelpts,nparams=1, nthetas=3){
  #browser()
  res <- .C("callEstimate",
            as.double((model$xmodel)),
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
  plot(f1$emulatedx, f1$emulatedy)
  print(f1)
  f2<-callEmulate(model, ans,m, nemupts=10)
  print(f2)
}

## use a given set of thetas to emulate the code
callEmulate <- function(model, thetas, nmodelpts, nparams=1, nthetas=3, nemupts=20, rangemin=0.0, rangemax=1.0){
  #browser()
  res <- .C("callEmulate",
            as.double((model$xmodel)),
            as.integer(nparams),
            as.double(model$training),
            as.integer(nmodelpts),
            as.double(thetas),
            as.integer(nthetas),
            finalx = double(nemupts),
            as.integer(nemupts),
            finaly = double(nemupts),
            finalvar = double(nemupts),
            as.double(rangemin),
            as.double(rangemax))
            
   results <- list(emulatedx=res$finalx, emulatedy=res$finaly, emulatedvar=res$finalvar)           
  results
}


callEmulateAtPoint <- function(model, thetas, point, nmodelpts, nparams=1, nthetas=3){
  #browser()
  res <- .C("callEmulateAtPt",
            as.double((model$xmodel)),
            as.integer(nparams),
            as.double(point),
            as.double(model$training),
            as.integer(nmodelpts),
            as.double(thetas),
            as.integer(nthetas),
            finaly = double(1),
            finalvar = double(1)
            )
  results <- list(des=point, mean=res$finaly, var=res$finalvar)
}

callEmulateAtList <- function(model, thetas, pointList,nemupts, nmodelpoints, nparams=1, nthetas=3){
  #browser()

  ## if(nemupts != dim(pointList)[1]){
  ##   cat("nemupts ", nemupts, "\n")
  ##   print("error not enough emulate points")
  ## }

  ## cat("calling Emulate at list\n")
  ## cat("nmodelpoints: ", nmodelpoints, "\n")
  ## cat("nemupts: ", nemupts, "\n")
  ## cat("nparams: ", nparams, "\n")
  
  res <- .C("callEmulateAtList",
            as.double((model$xmodel)),
            as.integer(nparams),
            as.double(pointList),
            as.integer(nemupts),
            as.double(model$training),
            as.integer(nmodelpoints),
            as.double(thetas),
            as.integer(nthetas),
            finaly = double(nemupts),
            finalvar = double(nemupts)
            )
  results <- list(des=pointList, mean=res$finaly, var=res$finalvar)
}


callEvalLikelyhood <- function(model, nmodelpoints, vertex, nparams=1,nthetas=4){
  answer <- 0.0
  likely <- .C("callEvalLikelyhood",
               as.double((model$xmodel)),
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
                                              


## the idea is to take a model in a few (ideally-inde) y-dimensions (factored into principle cpts or whatever in another fn)
## and emulate/estimate each dimension as if it was a single scalar field over the parameters.

# nydims -> the number of separate training dimensions
# training -> matrix of the y-values for each separate dimension
# model-> the parameter space over which we've evaluated all of the training vectors
# nmodelpts -> the length of each training vector, the number of points at which the model was evaluated
multidim <-  function(model, nmodelpts, training, nydims){
  nparams <- ncol(as.matrix(model))
  nthetas <- nparams+2 # this is the correct number for the gaussian cov fn
  bigthetas <- array(0, dim=c(nydims, nthetas)) # store all the thetas
  for(i in 1:nydims){ # estimate the thetas for each sub-model
    # we have to craft a custom frame for each call
    if(nydims > 1){
      bigthetas[i,] <-callEstimate(list(xmodel=model, training=training[,i]) 
                                   , nmodelpts, nparams, nthetas)
    } else {
      bigthetas[i,] <- callEstimate(list(xmodel=model, training) 
                                    , nmodelpts, nparams, nthetas)
    }
  }
  bigthetas
}


# passes the same set of training data twice, just to test that multidim makes some
# kind of sense
testMultiDim <- function(){
  source("~/Projects/Emulator/src/libRbind/testRbind.R")
  npts <- 10
  modelcomb <- demoModel(npts)
  thetas <- multidim(modelcomb$xmodel, npts, cbind(modelcomb$training, modelcomb$training), 2)
  thetas
}
