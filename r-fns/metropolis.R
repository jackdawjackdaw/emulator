# a test of some metropolis hastings integration

f  <- function(x){
  2*exp(-4*x*x)
}

q <- function(x, y){
  rnorm(1, x, y)
}

methaststep <- function(xt, rangemin, rangemax){
  yt <- runif(1, rangemin, rangemax) 
  prand <- runif(1, 0, 1)
  # the acceptance ratio
  paccept <- (f(yt)/f(xt))*(runif(1,rangemin, rangemax)/runif(1,rangemin, rangemax))
  # do the min
  if(paccept > 1.0){
    paccept <-1.0
  }
  if(prand <= paccept){
    xnext <- yt
  } else {
    xnext <- xt
  }
  xnext
}

loop <- function(nsteps){
  x0 <- runif(1, 0, 1)
  ans <- rep(NA, nsteps)
  for(i in 1:nsteps){
    x0 <- methaststep(x0, 0.0, 1.0)
    ans[i] = x0
  }
  ans
}

