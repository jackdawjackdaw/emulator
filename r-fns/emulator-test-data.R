library(lhs)

                                        
# generate a sample data set that can be used to test the emluator
# this generates a noisy flatline which snaps into a messy sin/cos wave
# at 35.7
test1d <- function(outname, numbersamples){
  answer <- rep(NA, numbersamples)
  x <- seq(0, 60, length=numbersamples)
  for(i in 1:numbersamples){
    answer[i] <- t1d(x[i]) + rnorm(1,0, 0.1)
  }
  write("# 1d test inputfile for emualtor", file=outname, append=FALSE)
final <- data.frame(x=x, y=answer)
write.table(final, file=outname, append=TRUE, sep=" ", col.names=FALSE, row.names=FALSE)
}

test1dRCall <- function(numbersamples){
  answer <- rep(NA, numbersamples)
  x <- seq(0, 60, length=numbersamples)
  for(i in 1:numbersamples){
    answer[i] <- t1d(x[i]) + rnorm(1,0, 0.1)
  }
  data.frame(xmodel=x, training=answer)
}


 t1d <- function(x){
   temp <- theta(x, 35.75)
   #print(temp)
   (sin(pi*x/5)+(1/5)*cos(4*pi*x/5))*temp
 }

 theta <- function(x, y){
   if( x >= y){
     1.0
   } else {
     0.0
   }
 }

# another test function, for 2d data
# y(x) = x1*exp(-x1^2 - x2^2)
# this is ranges on [-2,6], [-2, 6]
test2dExpo <- function(outname, numbersamples){
  samplepoints <- maximinLHS(numbersamples, 2)
  samplepoints <- (samplepoints * 8 - 2)
  temp <- rep(NA, numbersamples)
  for(i in 1:numbersamples){
    temp[i] <- y2d(samplepoints[i,]) + rnorm(1,0, 0.01)
  }
  answer <- data.frame(x=samplepoints[,1], y=samplepoints[,2], ans=temp)
  write("# 2d test inputfile for emulator", file=outname, append=FALSE)
  write.table(answer, file=outname,append=TRUE, sep=" ", col.names=FALSE, row.names=FALSE)
}

test2dExpoRCall <- function(numbersamples){
  samplepoints <- maximinLHS(numbersamples, 2)
  samplepoints <- (samplepoints * 8 - 2)
  temp <- rep(NA, numbersamples)
  for(i in 1:numbersamples){
    temp[i] <- y2d(samplepoints[i,]) + rnorm(1,0, 0.01)
  }
  answer <- data.frame(xmodel=samplepoints, training=temp)
}  

y2d <- function(x){
  x[1]*exp(-x[1]^2- x[2]^2)
}
