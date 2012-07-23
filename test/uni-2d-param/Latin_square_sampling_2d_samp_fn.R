#!/usr/bin/Rscript

## Run this program with either:
## 	$ R --slave < Latin_square_sampling_2d_samp_fn.R
## or:
##	$ chmod + x Latin_square_sampling_2d_samp_fn.R
## 	$ ./Latin_square_sampling_2d_samp_fn.R

## output file format
##	2 #d
##	200 #n
##	1 .. 400 #x
##	1 .. 200 #y
## The file format wil be perfect for input into interactive_emulator.


rangeMin  <- 0.0
rangeMax  <- 1.0
nmodelpts <- 200 ## NUMBER OF TRAINING POINTS
nparams   <- 2   ## number of dimensions
noutputs  <- 1   ## number of output scalars

multivariate_mode <- TRUE
#multivariate_mode <- FALSE

otype <- if (multivariate_mode) "multivar_" else "singlevar_"
file_name <- paste(otype, nmodelpts, "pts_",
	nparams, "d_", noutputs, "t.dat", sep="")

## set filename to "" for stdout.
file_name <- ""

if (! require("lhs", quietly=TRUE)) {
  R_LIBS_USER <- Sys.getenv("R_LIBS_USER")
  dir.create(R_LIBS_USER, FALSE, TRUE, "0777")
  REPO='http://cran.us.r-project.org'
  install.packages("lhs",lib=R_LIBS_USER,repos=REPO)
  cat('...try again...')
  q('no',1) # quit with return code 1
}

##
## this is an example function to be modeled.
##
modelFunc <- function(v) {
  x <- v[1];
  y <- v[2];
  (x-0.5)^2 + (y-0.5)^2
}

design <- (maximinLHS(nmodelpts, nparams) * (rangeMax - rangeMin)) + rangeMin
ymodel <- t(apply(design, 1, modelFunc))

fp <- if (nchar(file_name)) file(file_name,"w") else "";

options(digits=17) ## good enought to represent any double.

if (multivariate_mode)
  cat(noutputs, file=fp, fill=TRUE );
cat( nparams,   file=fp, fill=TRUE );
cat( nmodelpts, file=fp, fill=TRUE );
cat( t(design), file=fp, fill=TRUE );
cat( t(ymodel), file=fp, fill=TRUE );

if (nchar(file_name)) close(fp);
