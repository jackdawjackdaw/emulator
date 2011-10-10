## simple-regression example (1d)
## ccs, cec24@phy.duke.edu
##
## A simple example using R interface to libRBIND to make predictions of a function.
##
## We generate a design of n points in 1 dimension D = {u_1,..,u_n}
## using this we make n (slightly noisy) samples Y_t={Y_m(u_1),...Y_m(u_n)} of the trivial model yM(u) = 5*exp(-3*u)*(sin(u*10)) + 2.
##
## The samples and design are used to generate the optimal emulator hyper-parameters. The collection M={D, Y_T} is referred to 
## as the model in the code. The R interface requires that the model is a list with the design in a list-slot xmodel
## and the training evaluations in a list slot-training
##
## eg for some set of design points Design.points and realizations of Y_m Training.values
## model <- list(xmodel=Design.points, training=Training.values)
##
## For a power exponential covariance function C(u_i, u_j) = theta_0 * exp(-1/2 ( u_i - u_j) ^2 / theta_l^2) + theta_1.
## where in one parameter dimension, we have 
## 
## theta_0 overall scale
## theta_1 nugget
## theta_2 = theta_l sets the correlation length in the u space.
##
## To estimate the correct values for these hyper parameters we use the R function callEstimate.
##
## The set of estimated hyper parameters Theta, along with Y_t and D fully specify our emulator, we can now
## make predictions at a set of new points. We generate a set of points and pass this to callEmulateAtList, this returns 
## an list which contains the emulated mean and variance at each location.
##
## Finally we plot all this nicely.
##
## The numbered comments mark the logical flow of the above analysis.


## 1) load EmuRbind, this defines the R<->C interface, all the user callable emulator functions are here.
## 
source("~/local/include/libRbind/EmuRbind.R") # load the emu bindings

## 2) load libRBIND into R. Currently the codebase is NOT setup as an R package so this is a bit non standard.
arch <- system("uname -s", intern=TRUE)
if(is.loaded("callEstimate") == FALSE){
  libNAME <- "~/local/lib/libRBIND"
  if(arch == "Linux"){
    libNAME <- paste(libNAME, ".so", sep="")
  } else if(arch =="Darwin"){
    libNAME <- paste(libNAME, ".dylib", sep="")
  } else {
    buffer <- paste("error: uname -s gives ", arch , " not supported", sep="")
    stop(buffer)
  }
  dyn.load(libNAME)
}


## our simple model function
## with a little bit of noise
## the noise helps the estimator to converge in the limit of large N
yM <- function(z, delta=0.005) {5*exp(-3*z)*(sin(z*10)) + 2 + rnorm(1, mean=0, sd=delta) }

## construct the model data for emulation using a latin hypercube sampling method to generate the design
## @return a list containing the model design and the evaluations of the model over the design
make.model <- function(m, rangeMin=0.0, rangeMax=1.0, modelFunc=yM){
  retval <- require("lhs", quietly=TRUE)

  ## test for the library
  if(retval == FALSE){
    print("installing the lhs library")
    install.packages(c("lhs"))
  }
  ## we use the R library LHS to construct a latin hypercube design in 1d
  design <- (maximinLHS(m, 1))*(rangeMax-rangeMin) + rangeMin
  ymodel <-  modelFunc(design)  
  ## the model canonically contains an xmodel and a training set list, other entries
  ## are ignored
  model  <- list(xmodel = design, training = ymodel)
  invisible(model)
}


## estimate the hyperparameters for the supplied model
## @return a list containing the model, the estimated thetas and all the parameters needed to
## generate predictions from the model at some new untried input
estimate.model <- function(model){
  nmodelpts <- length(model$xmodel)   ## number of design points
  nparams <- 1 ## number of dimensions in the u space
  ## for the default power-exponential covariance function this nthetas = nparams + 2
  ## callEstimate will complain if this is set to an incorrect value.
  nthetas <- nparams + 3
  ## pick which covariance function to use from (1=powerExponential, 2=matern3/2, 3=matern5/2), 
  ## the power exponential function should almost always be used.
  cov.fn <- 4
  ## pick the order of the prior regression process (0=constant only, 1=linear model, 2=quadratic model, 3=cubic)
  reg.order <- 1

  ## callEstimate uses libRBIND to try and generate the best set of length scales for the supplied model
  ##
  ## note that model must contain xmodel and training slots which contain the model design and the 
  ## evaluations of the model at these design points
  ##
  ## IMPORTANT: models with multidimensional parameter spaces xmodel should be a matrix with nparams columns, 
  ## each row of the matrix specifies a location in the nparams-d design space.
  ##
  ## @return thetas.est a vector of the estimated covariance scales
  thetas.est <- callEstimate(model, nmodelpts,
                             nparams=nparams, nthetas=nthetas,
                             cov.fn=cov.fn, reg.order=reg.order)

  ## return a list with the appropriate thetas, the model data used to generate them, and
  ## which cov.fn and reg.order we used. Making predictions using the same thetas, and model but with a different
  ## cov.fn or reg.order than that used in the estimation process WILL NOT WORK! At best you'll crash the code,
  ## at worst you'll get some very confusing output.
  ##
  estim.result <- list(thetas=thetas.est, model=model, cov.fn=cov.fn, reg.order=reg.order,
                       nparams=nparams, nmodelpts=nmodelpts, nthetas=nthetas)
  invisible(estim.result)
}

## generate predictions for the model using the results of estimate.model
## @param estim.result list returned from estimate.model
## @return a list the points we evaluated the model at and the predicted mean and variance at each of these
predict.model <- function(estim.result){
  nemupts <- 128

  ## generate a uniform list of points across our 1d parameter space
  xmax <- max(estim.result$model$xmodel)
  xmin <- min(estim.result$model$xmodel)
  pointList <- seq(from=xmin, to=xmax, length.out=nemupts)

  ## callEmulateAtList uses libRBIND to generate predictions for the supplied model at
  ## the set of pointList locations in the parameter space.
  ## 
  ## we specify the model with the results from estimate.model
  emu.result <- callEmulateAtList(estim.result$model, estim.result$thetas,
                                  pointList, nemupts, nmodelpoints=estim.result$nmodelpts,
                                  nparams=estim.result$nparams, nthetas=estim.result$nthetas,
                                  reg.order=estim.result$reg.order,
                                  cov.fn=estim.result$cov.fn)
  invisible(emu.result)
}

## plot the model and the emulator results
## @param emu.result the list of emulator predictions generated by predict.model
## @param model the initial model created by make.model
plot.model <- function(emu.result, model){
  plot(model$xmodel, model$training, col="green", pch=2, ylim=c(0,5.5),
       xlab="u", ylab="Y")
  lines(emu.result$des, yM(emu.result$des), col="black",lty=2)
  lines(emu.result$des, emu.result$mean, col="red") ## plot the emulated mean and 95% confidence intervals
  lines(emu.result$des, emu.result$mean+1.96*sqrt(emu.result$var), col="red", lty=2)
  lines(emu.result$des, emu.result$mean-1.96*sqrt(emu.result$var), col="red", lty=2)

  legend("topright", c("true function", "training points", "emulated mean", "95% confidence bounds"),
         pch=c(-1,2, -1, -1), col=c("black", "green", "red", "red"), lty=c(2,-1, 1, 2))
}


## 3) generate the model with 10 samples of our test function
## re-run the script with different values here to see how changing the number of samples of our underlying model changes
## the shape of  our predictions
model <- make.model(10)

## 4) train the emulator by estimating the appropriate characteristic length scales (hyper parameters vector Theta)
estimResult <- estimate.model(model)
## estimResult now contains everything we need to make predictions, the estimation process is often slow and
## gives slightly different (but similar) results each time, if we want to make a set of consistent predictions now and
## at some later date it makes sense to save this object using R's "save" command and then we can use it again later
save(estimResult, file="estim-result-simple-example.dat")

## 5) make predictions across the interval spanned by our training set
emuResult <- predict.model(estimResult)

## 6) plot our predictions
plot.model(emuResult, model)

## this is hopefully enough to make predictions from a 1d model, for a higher dimensional parameter space the only
## things that should change are:
## xmodel needs to become a matrix with each row being a point in the design space
## nparams needs to be set to the dimension of the parameter space
## nthetas may need to be adjusted
## the pointList in predict.model needs to be expanded to multiple dimensions the R function expandGrid is very nice for this
## you'll need to figure out a way to plot a n-dim model (plot.model won't gracefully expand), 
## the R image and contour plot functions are pretty good.
## 
## enjoy...
