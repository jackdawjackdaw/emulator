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
write.table(answer, file=outname, append=TRUE, sep=" ", col.names=FALSE, row.names=FALSE)
}

 t1d <- function(x){
   temp <- theta(x, 35.75)
   print(temp)
   (sin(pi*x/5)+(1/5)*cos(4*pi*x/5))*temp
 }

 theta <- function(x, y){
   if( x >= y){
     1.0
   } else {
     0.0
   }
 }
