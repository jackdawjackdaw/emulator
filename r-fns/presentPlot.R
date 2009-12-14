## all the plots needed for the presentation
source("EmuRbind.R")
source("testRbind.R")

## run this to build all the plots
printPlots <- function(){
  postscript("demoModel-pts-only.ps")
  plotPlainModel()
  dev.off()
  postscript("demoModel-emulated.ps")
  plotEmulatedModel()
  dev.off()
  postscript("demoModel-emu-interp.ps")
  demoInterpolation()
  dev.off()
}

plotPlainModel <- function(){
  ## make the model data
  model <- demoModel(8, rangeMin=0.0, rangeMax=1.5)
  setEmulatorOptions(0,1.9)
  nemupts <- 200

  ## doesn't do any emulation
  sequence <- seq(0.0,1.5, length=nemupts)
  actual <- data.frame(x=sequence, y=yM(sequence))

  plot(model$xmodel, model$training, ylim=range(0.0, 6.0), xlab="x", ylab="y", pch=19, xlim=range(0.0,1.5))
  lines(actual, type="l", lwd=2, lty=2, col="darkolivegreen")
  title(main="y=5*exp(-3x)sin(10x)+2")
  legend(x=1.1, y=5, legend=c('training points','model'),
         col=c('black', 'darkolivegreen'),
         lwd=2,
         lty=c(0,2))
}

plotEmulatedModel <- function(){
  ## make the model data
  nmodelpts <- 8
  model <- demoModel(nmodelpts, rangeMin=0.0, rangeMax=1.5)
  setEmulatorOptions(0,1.9)
  nemupts <- 200

  ## doesn't do any emulation
  sequence <- seq(0.0,1.5, length=nemupts)
  actual <- data.frame(x=sequence, y=yM(sequence))
  ans <- callcode(model, nmodelpts, nemupts=nemupts, rangemin=0.0, rangemax=1.5)
  
  plotResultsTest(model, ans, actual, "y=5*exp(-3x)sin(10x)+2")
  

  legend(x=1.1, y=5, legend=c('Training points','Model', 'Emulated mean', 'Confidence'),
         col=c('black', 'darkolivegreen', 'red', 'red'),
         lwd=2,
         lty=c(0,2,1,2),
         pch=c(19,-1,-1,-1))
  
}  

