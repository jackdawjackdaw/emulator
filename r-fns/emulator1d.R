# the emulator in 1d, as a function (lots of copying crap around i know...)
library("mlegp")

#you need to set a seed before you do any of the fitting

#
# expects the input file to be in the form
# param1...paramN Y(params)
doEmulator <- function(inputFilename, npoints,rangeMin, rangeMax, outputFilename){
  # first we read the data in
  modeldata <- read.table(inputFilename, head=FALSE, comment.char="#")
  m <- nrow(modeldata) # no of points where the model was evaluated
  dim <- ncol(modeldata) -1
                                         
  params <- cbind(modeldata)[,1:dim] # every col but the last one
  yMdM <- cbind(modeldata)[,(dim+1)] # the penultimate col only
  #ySd <- cbind(modeldata)[,(dim+2)] # the last column

  ## maximum likelyhood
  ## introduce a nugget which is the variance of the data points
  mleYm <- mlegp(params, yMdM, constantMean=0)
  sig2M <- mleYm$sig2 # variance
  muM <- mleYm$mu # mean

  ## covariance
  ## the correlation fn here
  cM <- function(z1, z2){
    exp(-sum(mleYm$beta*(z1-z2)^2))
  }

  # a vector of the correlation of a given point
  # against the model points
  rz <- function(z){
    temp <- rep(-1, length=m)
    for(i in 1:m) temp[i] <- cM(z, params[i])
    return(sig2M*temp)
  }

                                        # make the covariance array
  covM <- diag(m) # empty
  for(j in 1:(m-1))
    for(i in (1+j):m)
      covM[i,j] <- covM[j,i] <- exp( - sum(mleYm$beta * (params[j] - params[i])^2))

### emulator fns for one input
  # for inputs z \nelemntof D^M, the mean can be estimated by: (cf. 4.4 & following mle estimate discussion)
  predMeanMLE <- function(z) { # z={t, u}
    rzlocal <- rz(z)
                                        # this is weird, why is yMdM being like this?
    return(muM[1] + as.numeric(rzlocal %*% solve(sig2M*covM, yMdM - muM)))
  }

  # measuring the uncertainty in the above approximation: (cf. 4.5f)
  # since the uncertainty of sig2M and mleYm$beta is not taken into account, this is an underestimate!
  predVarMLE <- function(z) { # this is the variance at the input z
    rzl <-rz(z)
    return (sig2M - as.numeric(rzl %*% solve(sig2M*covM, rzl)) )
  } # comment: this is supprisingly small... is it correct?

  predCovMLE <- function(z1,z2) { # this is the correlation between 2 random inputs
    rzl1 <-rz(z1)
    rzl2 <-rz(z2)
    return(sig2M*cM(z1, z2) - as.numeric(rzl1 %*%  solve(sig2M*covM, rzl2)))  
  }

  # now actually calculate the stuff
  x <- seq(rangeMin, rangeMax, length=npoints)

  emuMean <- rep(0, npoints)
  emuVar <- rep(0, npoints)
  emuVarSelf <- rep(0, npoints)

  for(i in 1:npoints){
    emuMean[i] <- predMeanMLE(x[i])
    emuVar[i] <- predVarMLE(x[i])
  }

  # create a results matrix
  results <- matrix(nrow = npoints, ncol = 3)
  for(i in 1:npoints){
    results[i, ] <- c(x[i], emuMean[i], emuVar[i])
  }

  # output to the correct file
  colnames <- c("#temp emuMean emuVar")
  useful <- sprintf("#infile = %s\n", inputFilename)
  write(useful, file=outputFilename, append=FALSE)
  write(colnames, file=outputFilename, sep=" ", append=TRUE)
  write.table(results, file=outputFilename, col.names =FALSE, row.names =FALSE, append=TRUE)
}
  
  


  
