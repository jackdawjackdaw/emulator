# a test of some metropolis hastings integration

f  <- function(x){
  exp(-4*x*x)
  #sin(3*x)
}

q <- function(x, y){
  rnorm(1, x, y)
}

# so this is an U[0,1] independance sampler. why doesn't it work?
mhsUniform <- function(xt, rangemin, rangemax){
  yt <- runif(1, rangemin, rangemax) 
  prand <- runif(1, 0, 1)
  # the acceptance ratio
  paccept <- (f(yt)/f(xt))*(runif(1,rangemin, rangemax)/runif(1,rangemin, rangemax))
  # do the min
  if(paccept > 1.0){
    #print("topout")
    paccept <-1.0
  }
  if(prand <= paccept){
    xnext <- yt
  } else {
    xnext <- xt
  }
  xnext
}

qUniWalk<- function(x, sigma){
  q <- x + sigma*runif(1, -1, 1) 
  q
}

mhsUniformWalk <- function(xt,sigma){
  #yt <- xt + sigma*runif(1, -1, 1)
  yt <- qUniWalk(xt, sigma)
  #print(yt)
  prand <-runif(1,0,1)  
  #paccept <- (f(yt)*qUniWalk(xt,sigma))/(f(xt)*qUniWalk(yt,sigma))
  paccept <- (f(yt)/f(xt))
  #print(paccept)
  if(paccept > 1.0){
    #print("topout")
    paccept <- 1.0
  }
  if(prand <= paccept){
    xnext <- yt
  }
  else {
    xnext <- xt
  }
  xnext
}
    
loop <- function(nsteps){
  ans <- rep(NA, nsteps)
  x <- runif(1, 0, 1)
  for(i in 1:nsteps){
    x<- mhsUniformWalk(x, 0.05)
    #x <- mhsUniform(x, -1.0, 1.0)
    ans[i] <- x
  }
  ans
}

test <- function(){
  nsamples <- 100000
  samples <- loop(nsamples)
  t <- hist(samples,40)
  plot(t$mids, (t$counts/max(t$counts)))
  #lines(density(samples))
  lines(t$mids, f(t$mids))
}
  

montegrate <- function(nsteps){
  ans <- rep(NA, nsteps)
  dif <- rep(NA, nsteps)
  x <- runif(1, 0, 1)
  
  for(i in 1:nsteps){
    ans[i] <- f(x)
    diff[i] <- abs(0.663331 - ans[i])
    x<- runif(1,0,1)
  }

  answer <- data.frame(samples=ans, diff=diff)
  answer
}
