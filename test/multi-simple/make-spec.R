nparams <- 3
noutputs <- 6

design <- read.table("./data/design_2_sorted.dat")
lumfn <- read.table("./data/lum_fun_outputs_2.dat")
metallicity <- read.table("./data/metallicity_MV_outputs_2.dat")


outputs <- cbind(lumfn, metallicity)

nmodelpoints <- dim(outputs)[1]

## now write the stuff to stdout
cat(noutputs, "\n")
cat(nparams,  "\n")
cat(nmodelpoints, "\n")

for(i in 1:nmodelpoints){
  for(j in 1:nparams){
    cat(design[i,j], " ")
  }
  cat("\n")
}

for(i in 1:nmodelpoints){
  for(j in 1:noutputs){
    cat(outputs[i,j], " ")
  }
  cat("\n")
}
