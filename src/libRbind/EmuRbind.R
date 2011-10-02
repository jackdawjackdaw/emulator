library("lhs")

## ccs, cec24@phy.duke.edu (2011-2008)
##
## functions for calling the emulator through R
##
## common parameters
## 
## @param model: a list/data frame with entries: xmodel & training
## @param model$xmodel: the design positions where the model was evaluated
## (nparams cols, nmodelpoints rows)
## @param model$training: the values of the model at each of the positions in xmodel
## @param nparams: number of design dimensions
## @param nthetas: number of hyperparams for the emulator (needs to agree with the cov-fn) 
##
## @param reg.order sets the order of the linear regression model applied,
## where 0: a0, 1: a0 + a1*x etc
##
## @param cov.fn {1,2,3} selects the covariance function from { power_exp, matern_32, matern_52 }
## the cov.fns require nthetas as such:
## power_exp -> nparams + 2 
## matern_32 -> 3
## matern_52 -> 3
## 
## the thetas must always be stored as: {amplitude, nugget, length_scales...}
## amplitude -> overall scale applied to the cov fn
## nugget -> (iid noise term) applied to the diagonal
## length_scales -> for power_exp there are nparams of these one for each direction,
## for the matern class there is a single length scale
##
##
## note: we can somewhat debug the running code by starting r with gdb underneath it as
## R -d gdb --vanilla
## this will launch gdb, you can then run R from gdb, load the source file containing the
## offensive code and when things go screwy you'll have gdb underneath to give some clues
## 

 
## estimates the thetas for a model, do this first to create your emulator
## 
## @param fixedNugget: force the nugget to be within 5% of the given value, don't if NULL
## @return "optimal" estimated thetas for the supplied model, cov.fn, reg.order combo
##
## the created emulator is really the set: {model, thetas (which come from here) and reg.oder and cov.fn)
## so to use the emulator to make predictions at locations you'll need to send the same {model, thetas, reg.order and cov.fn)
## 
## 
callEstimate <- function(model, nmodelpts,nparams=1, nthetas=3, fixedNugget=NULL, cov.fn=1, reg.order=1){
  checkCovFn(nthetas, nparams, cov.fn)

  if(is.null(fixedNugget)){
    res <- .C("callEstimate",
              as.double(model$xmodel),
              as.integer(nparams),
              as.double(model$training),
              as.integer(nmodelpts),
              as.integer(nthetas),
              thetas = double(nthetas),
              as.integer(0), as.double(0),
              as.integer(cov.fn), as.integer(reg.order))
  } else {
    res <- .C("callEstimate",
              as.double(model$xmodel),
              as.integer(nparams),
              as.double(model$training),
              as.integer(nmodelpts),
              as.integer(nthetas),
              thetas = double(nthetas),
              as.integer(1), as.double(fixedNugget),
              as.integer(cov.fn), as.integer(reg.order))
  }
  res$thetas
}

## 
## the idea is to take a model in a few (ideally-inde) y-dimensions (factored into principle cpts or whatever in another fn)
## and emulate/estimate each dimension as if it was a single scalar field over the parameters.
# @param nydims -> the number of separate training dimensions
# @param training -> matrix of the y-values for each separate dimension
#
# @param fixedNugget see above
#
# @return matrix of optimal thetas one row per training dimension
#
multidim <-  function(model, nmodelpts, training, nydims, fixedNugget=NULL, cov.fn=1, reg.order=1){
  nparams <- ncol(as.matrix(model))
  
  if(cov.fn == 1){
    nthetas <- nparams+2 # this is the correct number for the gaussian cov fn
  } else {
    nthetas <- 3
  }

  ## store all the thetas and pass them one row at a time
  bigthetas <- array(0, dim=c(nydims, nthetas)) 

  checkCovFn(nthetas, nparams, cov.fn)

  for(i in 1:nydims){ # estimate the thetas for each sub-model
    # we have to craft a custom frame for each call
    if(nydims > 1){
      bigthetas[i,] <-callEstimate(list(xmodel=model, training=training[,i]) 
                                   , nmodelpts, nparams, nthetas, fixedNugget[i],
                                   cov.fn=cov.fn, reg.order=reg.order)
    } else {
      bigthetas[i,] <- callEstimate(list(xmodel=model, training) 
                                    , nmodelpts, nparams, nthetas, fixedNugget[i],
                                    cov.fn=cov.fn, reg.oder=reg.order)
    }
  }
  bigthetas
}




## return the emulated mean and var at point
## N calls of this are substantially slower than using callEmulateAtList,
## the AtList function is much more efficient.
#
# @param point the locn to emulate the model at
# 
# @return a list res containing fields {des, mean, var}
# @return res$des location at which we obtain the emulated mean and variance
# @return res$mean emulated mean at res$des
# @return res$var emulated variance at res$des
callEmulateAtPoint <- function(model, thetas, point, nmodelpts, nparams=1, nthetas=3,
                               cov.fn=1, reg.order=1){

  checkCovFn(nthetas, nparams, cov.fn)
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
            finalvar = double(1),
            as.integer(cov.fn),
            as.integer(reg.order)
            )
  results <- list(des=point, mean=res$finaly, var=res$finalvar)
}


## the preferred way to compue the emulated mean and variance at a set of points given by pointList
##
## @param pointList an R matrix (nrow=nemupts, ncol=nparams) of points at which to emulate the model
## @param nemupts how many points to compute
## @return a list res(des, mean, var)
## @return res$des<-pointList
## @return res$mean emulated mean at each location in pointList
## @return res$var emulated var at each location in pointList
callEmulateAtList <- function(model, thetas, pointList,nemupts, nmodelpoints, nparams=1, nthetas=3,
                              cov.fn=1, reg.order=1){

  checkCovFn(nthetas, nparams, cov.fn)
  
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
            finalvar = double(nemupts),
            as.integer(cov.fn),
            as.integer(reg.order)
            )

  results <- list(des=pointList, mean=res$finaly, var=res$finalvar)
}

## compute the lhood of a set of points in hyper-param space, given a set of thetas, a cov.fn and a regression
## model
## 
## @param pointList an R matrix (nrow=nevalPoints, ncol=nthetas) of the locns in the hyperparam space to
## estimate the log lhood
## @param nevalPoints how many locations
## @return a list res(des, lhood)
## @return res$des <- pointList
## @return res$lhood vector of lhoods at each location in pointList
callEvalLhoodList <- function(model, pointList, nevalPoints, nmodelPoints, nparams=1, nthetas=3,
                              cov.fn=1, reg.order=1){
  checkCovFn(nthetas, nparams, cov.fn)
  
  answer <- 0.0
  
  res <- .C("callEvalLhoodList",
            as.double(model$xmodel),
            as.integer(nparams),
            as.double(pointList),
            as.integer(nevalPoints),
            as.double(model$training),
            as.integer(nmodelPoints),
            as.integer(nthetas),
            finalLhood = double(nevalPoints),
            as.integer(cov.fn),
            as.integer(reg.order)
            )

  results <- list(des=pointList, lhood=res$finalLhood)
}

##
## not to be externally called, checks that the cov.fn and nthetas given match up
## should be called before doing any .C calls into rbind
checkCovFn <- function(nthetas, nparams, cov.fn){
  if(cov.fn < 1 || cov.fn > 3){
    buffer <- paste("cov.fn index out of range: ", cov.fn, sep="")
    stop(buffer)
  }

  if(cov.fn == 1){ ## power exp
    if(nthetas != nparams + 2){
      buffer <- paste("power exp called with nthetas: ", nthetas,
                      "needs nthetas: ", nparams+2, sep="")
      stop(buffer)
    }
  } else { ## matern
    if(nthetas != 3){
      buffer <- paste("matern class need nthetas = 3, supplied: ", nthetas)
      stop(buffer)
    }
  }
}

##
## fns for MC calls to the emulator

## first setup the structure in rbind.h with copies of the model, options and
## inverse cov matrix etc, so we don't have to do the inversion for each
## evaluation in the MC chain.
setup.MC <- function(model, thetas, nmodelpoints, nparams=1, nthetas=3, cov.fn=1, reg.order=1){
  checkCovFn(nthetas, nparams, cov.fn)

  res <- .C("setupEmulateMC",
            as.double((model$xmodel)),
            as.integer(nparams),
            as.double(model$training),
            as.integer(nmodelpoints),
            as.double(thetas),
            as.integer(nthetas),
            as.integer(cov.fn),
            as.integer(reg.order)
            )
}

## make a quick call to the emulator mean and variance given a previous setup
##
emulate.MC <- function(point){
  res <- .C("callEmulateMC",
            as.double(point),
            finalMean = double(1),
            finalVar = double(1))

  results <- list(des=point, mean=res$finalMean, var=res$finalVar)  
}

## clear the mc data, have to do this before making another call
destroy.MC <- function(){
  .C("freeEmulateMC")
}


##
## this sets up the MC for a multi-variate emulator, 
##
## follow the same procedure as multi-dim
setup.MC.multi <- function(model, thetas, nmodelpoints, nydims, nparams, nthetas, cov.fn=1, reg.order=1){
  checkCovFn(nthetas, nparams, cov.fn)

  res <- .C("setupEmulateMCMulti",
            as.double((model$xmodel)), ## matrix (nrows=nmodelpts, ncols=nparams)
            as.integer(nparams),
            as.double(model$training), ## matrix (nrows=nmodelpts, ncols=nydims)
            as.integer(nydims),
            as.integer(nmodelpoints),
            as.double(thetas), ## matrix (nrows=nydims, ncols=nthetas)
            as.integer(nthetas),
            as.integer(cov.fn),
            as.integer(reg.order))

}

## make a quick call to the emulator mean and variance given a previous setup
##
emulate.MC.multi <- function(point, nydims){
  res <- .C("callEmulateMCMulti",
            as.double(point),
            as.integer(nydims),
            finalMean = double(nydims),
            finalVar = double(nydims))

  results <- list(des=point, mean=res$finalMean, var=res$finalVar)  
}


## delete everything we setup
destroy.MC.multi <- function(nydims){
  .C("freeEmulateMCMulti",
     as.integer(nydims))
}
                         
  



