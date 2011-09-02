library("lhs")

## ccs, cec24@phy.duke.edu
##
## functions for creating and sampling the mean and var of an emulator over a data set
##
## 
##
## note: we can somewhat debug the running code by starting r with gdb underneath it as
## R -d gdb --vanilla
## this will launch gdb, you can then run R from gdb, load the source file containing the
## offensive code and when things go screwy you'll have gdb underneath to give some clues
## 

 
## estimates the thetas for a model 
# if fixedNugget is not set to NULL the supplied value is used to fix the nugget
# for the estimation process. This can be used to force some uncertainty in the training points
# also it can just fuck things up majorly
callEstimate <- function(model, nmodelpts,nparams=1, nthetas=3, fixedNugget=NULL){
  #browser()
  
  if(is.null(fixedNugget)){
    res <- .C("callEstimate",
              as.double(model$xmodel),
              as.integer(nparams),
              as.double(model$training),
              as.integer(nmodelpts),
              as.integer(nthetas),
              thetas = double(nthetas),
              as.integer(0), as.double(0))
  } else {
    res <- .C("callEstimate",
              as.double(model$xmodel),
              as.integer(nparams),
              as.double(model$training),
              as.integer(nmodelpts),
              as.integer(nthetas),
              thetas = double(nthetas),
              as.integer(1), as.double(fixedNugget))
  }
  res$thetas
}

## 
## the idea is to take a model in a few (ideally-inde) y-dimensions (factored into principle cpts or whatever in another fn)
## and emulate/estimate each dimension as if it was a single scalar field over the parameters.
# nydims -> the number of separate training dimensions
# training -> matrix of the y-values for each separate dimension
# model-> the parameter space over which we've evaluated all of the training vectors
# nmodelpts -> the length of each training vector, the number of points at which the model was evaluated
multidim <-  function(model, nmodelpts, training, nydims, fixedNugget=NULL){
  nparams <- ncol(as.matrix(model))
  nthetas <- nparams+2 # this is the correct number for the gaussian cov fn
  bigthetas <- array(0, dim=c(nydims, nthetas)) # store all the thetas

  
  for(i in 1:nydims){ # estimate the thetas for each sub-model
    # we have to craft a custom frame for each call
    if(nydims > 1){
      bigthetas[i,] <-callEstimate(list(xmodel=model, training=training[,i]) 
                                   , nmodelpts, nparams, nthetas, fixedNugget[i])
    } else {
      bigthetas[i,] <- callEstimate(list(xmodel=model, training) 
                                    , nmodelpts, nparams, nthetas, fixedNugget[i])
    }
  }
  bigthetas
}




## return the emulated mean and var at point
## N calls of this are substantially slower than using callEmulateAtList,
## the AtList function is much more efficient.
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


## the preferred way to compue the emulated mean and variance at a set of points given by pointList
callEmulateAtList <- function(model, thetas, pointList,nemupts, nmodelpoints, nparams=1, nthetas=3){
  
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

## compute the lhood of a set of points in hyper-param space
callEvalLhoodList <- function(model, pointList, nevalPoints, nmodelPoints, nparams=1, nthetas=3){
  answer <- 0.0
  
  res <- .C("callEvalLhoodList",
            as.double(model$xmodel),
            as.integer(nparams),
            as.double(pointList),
            as.integer(nevalPoints),
            as.double(model$training),
            as.integer(nmodelPoints),
            as.integer(nthetas),
            finalLhood = double(nevalPoints)
            )

  results <- list(des=pointList, lhood=res$finalLhood)
}

## deprecated fns here
callEmulate <- function(model, thetas, nmodelpts, nparams=1, nthetas=3, nemupts=20, rangemin=0.0, rangemax=1.0){
  stop("deprecated, use callEmulateAtPoint or callEmulateAtList")
}

callcode <- function(model, nmodelpts, nparams=1, nthetas=3, nemupts=50, rangemin=0.0, rangemax=4.0){
  stop("callcode is deprecated, use callEstimate and callEmulateAtList")
}
  
callEvalLikelyhood <- function(model, nmodelpoints, vertex, nparams=1,nthetas=3){
  stop("callEvalLikelyhood is deprecated, use callEvalLhoodList")
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

