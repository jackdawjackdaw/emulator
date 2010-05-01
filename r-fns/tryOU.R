## try and sample the OU process
## and then model that sucker.
source("EmuRbind.R")
source("plotCovReals.R")

#m <- 8

nbig <- 150

model <- makeReal(1.0, nbig, thetaLength=0.001)


emulateOU <- function(ouModel, m, color){

 design <- data.frame(xmodel=rep(NA, m), training=rep(NA, m))

 step <- floor(length(ouModel$x) / m )

 for(i in 1:m){
   design$xmodel[i] <- model$x[i*step]
   design$training[i] <- model$y[i*step]
 }

 setEmulatorOptions(0,1.9)
 #browser()
 ans <- callcode(design, m, nemupts=100, rangemin=0.0, rangemax=1.0)

 #legend(x=0.8,y=2.0, legend=c('actual process', 'training points'), col=c('black', 'red'), bg='white')
 points(design, pch=12, col=color)
 lines(ans$emulatedx, ans$emulatedy, lwd=2, col=color)
 confidence <- rep(NA, length(ans$emulatedvar))
 for(i in 1:length(ans$emulatedvar))
   confidence[i] <- sqrt(ans$emulatedvar[i])*1.65585

 lines(ans$emulatedx, ans$emulatedy + confidence, col=color, lty=2)
 lines(ans$emulatedx, ans$emulatedy - confidence, col=color, lty=2)

}


plot(model, ylim=range(floor(min(model$y))-0.25,floor(model$y)+0.5), pch="+")
title(main="Sampling the OU process, sigma=0.001")
lines(model, lwd=1)
grid()
emulateOU(model, 12, 'red')
#emulateOU(model, 15, 'darkolivegreen')

legend(x="bottomright", legend=c('model', 'emulator', 'confidence intervals'),
         col=c('black', 'red', 'red'),
         lwd=2,
         #pch=c(12, "+", -1, -1),
         lty=c(1,1,2),
         bg='white')
