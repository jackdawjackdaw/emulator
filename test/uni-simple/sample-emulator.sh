#!/bin/zsh
## ccs, generate samples from the trained model
##
## binpath reaches into build, is this a bad idea?
#binpath=../../build/src/interactive_emulator
binpath=~/local/bin/interactive_emulator
$binpath interactive_mode univariate_snapshot_file < sample_locations.dat > temp.dat
## now strip out the even/odd lines
awk '(NR % 2 == 0 && NR > 5)' temp.dat > emu_mean.dat
awk '(NR % 2  && NR > 5)' temp.dat > emu_var.dat

## clear up
if [ -f temp.dat ]; then
		rm temp.dat
fi

## the R script plot-results.R will then make a nice plot if you want
