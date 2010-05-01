source("EmuRbind.R")
source("testRbind.R")

m <- 5
#yM <- function(x) { 0.1*x^2*sin(5*x)}
                 

makeModel <- function(npts, rangeMin, rangeMax){
  x <- seq(from=rangeMin, to=rangeMax, length.out=npts)
  y <- yM(x)
  answer <- data.frame(xmodel=x, training=y)
  answer
}

joinModels <- function(modelLeft, modelRight){  
  length <- length(modelLeft$xmodel) + length(modelRight$xmodel)

  x <- append(modelLeft$xmodel, modelRight$xmodel)
  y <- append(modelLeft$training, modelRight$training)

  answer <- data.frame(xmodel=x, training=y)
  answer
}    
         



actualX <- seq(from=0.0, to=3.5, length.out=100)
actualY <- rep(NA, 100)
actualY <- yM(actualX)

model <- makeModel(15, 0.0, 3.5)

modelLeft <- makeModel(7, 0.0, 0.5)
modelRight <- makeModel(7, 1.0, 2.0)

modelJoined <- joinModels(modelLeft, modelRight)

plot(modelJoined, col="red", ylim=range(0.0, 6.0), xlim=range(0.0,2.0), pch=12)
lines(actualX, actualY, col='darkolivegreen')

setEmulatorOptions(0, 1.9)
ans <- callcode(modelJoined, 14, nemupts=100, rangemin=0.0, rangemax=3.5)


color="red"

lines(ans$emulatedx, ans$emulatedy, lwd=2, col=color)
 confidence <- rep(NA, length(ans$emulatedvar))
 for(i in 1:length(ans$emulatedvar))
   confidence[i] <- sqrt(ans$emulatedvar[i])*1.65585

 lines(ans$emulatedx, ans$emulatedy + confidence, col=color, lty=2)
 lines(ans$emulatedx, ans$emulatedy - confidence, col=color, lty=2)

legend(x="bottomright", legend=c('model', 'emulator', 'confidence intervals'),
         col=c('black', 'red', 'red'),
         lwd=2,
         #pch=c(12, "+", -1, -1),
         lty=c(1,1,2),
         bg='white')
