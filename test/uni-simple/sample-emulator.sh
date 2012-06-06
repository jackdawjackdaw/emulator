#!/bin/sh
## ccs, generate samples from the trained model
##
## binpath reaches into build, is this a bad idea?
binpath=../../build/src/interactive_emulator
$binpath interactive_mode univariate_snapshot_file 
