# test the linear interpolation model,
# it works ok but this is not a good demo yet
testLinInt <- function(m){
  model <- demoModel(m, 1)
  f <- lm(model$training ~ model$xmodel + I(model$xmodel^2) + I(model$xmodel^3) +I(model$xmodel^4))
  f1 <- lm(model$training ~ model$xmodel)
  plot(model, ylim=range(0,6))
  #print(summary(f) )
  #abline(f)
  x <- seq(0,1,by=0.01)
  y<-rep(0, 101)
  y<-f[1]$coefficients[1]*x^0 + f[1]$coefficients[2]*x + f[1]$coefficients[3]*(x^2) + f[1]$coefficients[4]*x^3 +f[1]$coefficients[4]*x^4
  #y1 <- rep(0, 101)
  #y1 <- f1[1]$coefficients[1]*x^0 + f1[1]$coefficents[2]*x

  lines(x,y)
  abline(f1, col="blue")
  ## ModelPlus <- demoModel(m+2, 1)
  ## points(ModelPlus, ylim=range(0,6), col="red")

  ## now run the emulator and save the day!
  setEmulatorOptions(0,1.9) ## set default ops gauss cov fn
  results <- callcode(model, m, rangemin=0.0, rangemax=1.5, nemupts=300)
  
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- sqrt(results$emulatedvar[i])*1.65585
  lines(results$emulatedx, results$emulatedy, col="red", lwd=2)
  
  title(main="Comparing Least Squares Regression")
  lines(results$emulatedx, results$emulatedy + confidence, col="red", lty=2)
  lines(results$emulatedx, results$emulatedy - confidence, col="red", lty=2)


}

cols=c("blue", "violet", "darkred", "deepskyblue")

demoInterpolation <- function(){
  x<- seq(0.0, 1.5, length.out=100)
  y<- yM(x)
  m <- 8
  master <- demoModel(m, 0, rangeMin=0.0, rangeMax=1.5)

  plot(master$xmodel, master$training, pch=19, col="black", xlim=range(0,1.5), ylim=range(0,6), cex=1.5, xlab="x", ylab="y")  
  ## now run the emulator and save the day!
  setEmulatorOptions(0,1.9) ## set default ops gauss cov fn
  results <- callcode(master, m, rangemin=0.0, rangemax=1.5)
  
  confidence <- rep(NA, length(results$emulatedvar))
  for(i in 1:length(results$emulatedvar))
    confidence[i] <- sqrt(results$emulatedvar[i])*1.65585

  ## colPoly <- rgb(190, 190, 190, alpha=70, maxColorValue=255)
  ## xxConf <- c(results$emulatedx, rev(results$emulatedx))
  ## yyConf <- c(results$emulatedy + confidence, rev(results$emulatedy - confidence))
  ## polygon(xxConf, yyConf, col=colPoly, lty=2)              

  # plot the emulated mean on top of the rest
  lines(results$emulatedx, results$emulatedy, col="red", lwd=2)
  
  title(main="Comparing polynomial interpolation against GP emulator")
  lines(results$emulatedx, results$emulatedy + confidence, col="red", lty=2)
  lines(results$emulatedx, results$emulatedy - confidence, col="red", lty=2)

  ## do the interp parts


  for(i in 1:3){
    order <- i*2+1
    ## just pull a little bit from the master
    ## this is not a very elegant way to do things...
    miniModel <- data.frame(xmodel=master$xmodel[1:order], training=master$training[1:order])
    #model <- demoModel(i, 0)
    #interp <- runLagInt(model)
    interp <- runLagInt(miniModel, interpPoints=200)
    lines(interp$x, interp$y, col=cols[i], lwd=2)

  }
  lines(x,y, col="darkolivegreen", lwd=3, lty=2)
  legend(x=1.0, y=5, legend=c('model', 'emulator',  'emulator 90% confidence',
                       '3rd order', '5th order', '7th order'),
         col=c('darkolivegreen', 'red', 'red', cols[1], cols[2], cols[3]),
         bg='white',
         lwd=2, cex=0.75,
         lty=c(2,1,2,1,1,1))
  grid()

}
    

## will be upset if you call this without an open plot
runLagInt <- function(model, rangeMin=0.0, rangeMax=1.5, interpPoints=100){
  interpX <- seq(rangeMin, rangeMax, length.out=interpPoints)
  interpY <- rep(NA, interpPoints)
  for(i in 1:interpPoints)
    interpY[i] <- callInterpolate(model$xmodel, model$training, interpX[i])
  result <- data.frame(x=interpX, y=interpY)
  result
}

