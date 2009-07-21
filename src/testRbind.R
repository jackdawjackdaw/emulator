dyn.load("libRbind.so")

nparams <- 1
nthetas <- 4
nmodelpts <- 10
nemupts <- 30
rangemin <- 0.0
rangemax <- 4.0

model <- read.table("../model.dat")



callcode <- function(){
  .C("callEmulator",
     as.double(model[,1]),
     as.integer(nparams),
     as.double(model[,2]),
     as.integer(nmodelpts),
     as.integer(nthetas),
     finalx = double(nemupts*nparams),
     as.integer(nemupts),
     finaly = double(nemupts),
     finalemu = double(nemupts),
     as.double(rangemin),
     as.double(rangemax))
}
   
   
