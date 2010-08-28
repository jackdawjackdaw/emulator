source("EmuRbind.R")
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
      bigthetas[i,] <-callEstimate(data.frame(xmodel=model, training=training[,i]) 
                                   , nmodelpts, nparams, nthetas)
    } else {
      bigthetas[i,] <- callEstimate(data.frame(xmodel=model, training) 
                                    , nmodelpts, nparams, nthetas)
    }
  }
  bigthetas
}


# passes the same set of training data twice, just to test that multidim makes some
# kind of sense
testMultiDim <- function(){
  source("testRbind.R")
  npts <- 10
  modelcomb <- demoModel(npts)
  thetas <- multidim(modelcomb$xmodel, npts, cbind(modelcomb$training, modelcomb$training), 2)
  thetas
}
  
