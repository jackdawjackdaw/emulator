# a test of some metropolis hastings integration

# cacuchy dist
f  <- function(x){
  #exp(-4*x*x)
  #sin(3*x)
  # the cauchy distribution
  (1.0)/(1+x*x)
}

# this is for integrating the function exp(-4*x^2)
g <- function(x){
  pi*(1+x^2)*exp(-4*x^2)*(x>0)*(x<1)
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

qUniWalk<- function(x, sigma, a=0, b=1){
  q <- x + sigma*runif(1, -1, 1)
  if(q > b || q < a){
    q <- a + (q - a)%%(b-a);
  }    
  q
}

mhsUniformWalk <- function(xt,sigma){
  #yt <- xt + sigma*runif(1, -1, 1)
  yt <- qUniWalk(xt, sigma, -Inf, Inf)
  #print(yt)
  prand <-runif(1,0,1)
  # symmetry
  #paccept <- (f(yt)*qUniWalk(xt,sigma))/(f(xt)*qUniWalk(yt,sigma))
  paccept <- (f(yt)/f(xt))
  if(prand <= paccept){
    xnext <- yt
  }
  else {
    xnext <- xt
  }
  xnext
}
    
loop <- function(nsteps, delta=0.05){
  ans <- rep(NA, nsteps)
  x <- runif(1, 0, 1)
  for(i in 1:nsteps){
    x<- mhsUniformWalk(x, delta)
    ans[i] <- x
  }
  ans
}

test <- function(){
  nsamples <- 100000
  samples <- loop(nsamples, 2.0)
  t <- hist(samples,40)
  plot(t$mids, (t$counts/max(t$counts)))

  #lines(density(samples))
  lines(t$mids, f(t$mids))
  samples
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
