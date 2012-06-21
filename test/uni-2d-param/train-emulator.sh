#!/bin/sh
## ccs, train an emulator on input_model_file.dat using a first order regression
##
## binpath reaches into build, is this a bad idea?
#binpath=../../build/src/interactive_emulator
binpath=interactive_emulator
inputfile=./Latin_square_sampling_2d_samp_fn_200.dat
outputfile=M.dat
## 
$binpath estimate_thetas $inputfile $outputfile --regression_order=1
