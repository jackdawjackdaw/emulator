#!/bin/zsh
## ccs, generate samples from the trained model
##
## binpath reaches into build, is this a bad idea?
binpath=../../build/src/interactive_emulator
modelfile=./M.dat
samplefile=./sample_locations.dat

nskip=6 ## 4 + nparams

$binpath interactive_mode $modelfile < $samplefile > temp.dat

## if this worked...
if [ -f temp.dat ]; then 
## now strip out the even/odd lines
## this is a bit annoying, if nskip is even we have to do this the other way around
		if [ $(( $nskip % 2 )) -eq 0 ]; then
				awk "(NR % 2 == 0 && NR > $nskip)" temp.dat > emu_var.dat
				awk "(NR % 2  && NR > $nskip)" temp.dat > emu_mean.dat
				
		else
				## odd 
				awk "(NR % 2 == 0 && NR > $nskip)" temp.dat > emu_mean.dat
				awk "(NR % 2  && NR > $nskip)" temp.dat > emu_var.dat
		fi

## and clear up
		#rm temp.dat
fi

## the R script plot-results.R will then make a nice plot if you want