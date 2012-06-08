#!/bin/sh

# 5 #d
# 100 #n
# 1 .. 500 #x
# 1 ..100 #y

RANGEMIN="0.0"
RANGEMAX="1.0"
NMODELPTS="200" # NUMBER OF TRAINING POINTS
NPARAMS="2"   # number of dimensions
FILE_NAME="Latin_square_sampling_${NPARAMS}d_samp_fn_${NMODELPTS}.dat"

R --slave --vanilla <<EOF
cat.matrix <- function(m, file="", width=FALSE) {
  for (i in 1:dim(m)[1]) {
    cat(m[i,], fill=width, file=file);
    if (! width)
      cat("\n", file=file);
  }
}
cat.vector <- function(v, file="", width=FALSE) {
  cat(v, fill=width, file=file);
  cat("\n", file=file);
}

rangeMin  <- $RANGEMIN
rangeMax  <- $RANGEMAX
nmodelpts <- $NMODELPTS
nparams   <- $NPARAMS

if (! require("lhs", quietly=TRUE)) {
  R_LIBS_USER <- Sys.getenv("R_LIBS_USER")
  dir.create(R_LIBS_USER, FALSE, TRUE, "0777")
  REPO='http://cran.us.r-project.org'
  install.packages("lhs",lib=R_LIBS_USER,repos=REPO)
  cat('...try again...')
  q('no',1) # quit with return code 1
}

design <- (maximinLHS(nmodelpts, nparams) * (rangeMax - rangeMin)) + rangeMin

#
# this is an example function to be modeled.
#
modelFunc <- function(v) {
  x <- v[1];
  y <- v[2];
  (x-0.5)^2 + (y-0.5)^2
}



#
#modelFunc <- function(v) {
#  x <- v[1]; y <- v[2];
#  x + 1.5 * y
#  100 * sin(20 * x) + 100 * sin(20*y + 10*x)
#  0.6 * exp( -32 * ((x - 0.5)^2 + (y - 0.5)^2)) + 0.2 * sin(20 * x) + 0.2 
#}
#modelFunc <- function(x) {
#  5 * exp(-3*sum(x*x)) * sin(x[3] * 10)
#  + 4 * sin(x[1] * 13 + x[2] * 7)  + 2
#  + rnorm(1, mean=0, sd=0.1)
#}
#modelFunc <- function(x) {
#  d = 0.5 - x;
#  return exp(-10 * sum(d*d)) + 0.25 * (sin(x[1] * 20) + sin(x[2] * 15));
#}

ymodel <- apply(design, 1, modelFunc)

#
# The file format wil be perfect for input into interactive_emulator.
#

fp <- file("${FILE_NAME}","w")

options(digits=17)
cat(nparams, file=fp);
cat("\n", file=fp);
cat(nmodelpts, file=fp);
cat("\n", file=fp);
cat.matrix(design, fp, 72)
cat.vector(ymodel, fp, 72)

close(fp)

EOF
test "$?" = 0 && echo "${FILE_NAME}" 