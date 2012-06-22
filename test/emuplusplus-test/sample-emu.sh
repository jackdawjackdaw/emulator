#!/bin/sh
## sample the emulator using our test program
if [ ! -e ./emu-dir/univariate_snapshot_file ]; then 
		echo "# no snapshot file, run train-emu.sh first"
		exit -1
fi

test_Emu++ ./emu-dir/univariate_snapshot_file
