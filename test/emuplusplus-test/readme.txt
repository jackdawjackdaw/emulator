22.06.2012
ccs, cec24@phy.duke.edu

A test / example of using the EmuPlusPlus lib

The EmuPlusPlus lib provides access to the class emulator which can be used to sample an already trained emulator from within a c++ program. After running build.sh the library and header should be installed to $INSTALL_PATH/lib and $INSTALL_PATH/include respectively. 

RUNNING/BUILDING 

1) train the emulator using the interactive_emulator code:
	 ./train-emu.sh
	 this creates a statefile in ./emu-dir called univariate_snapshot_file we will use this snapshot file to sample the emlator

2) build the example binary:
	 ./build.sh
	 this should compile and install test_Emu++, most likely to ~/local/bin

3) sample our emulator:
	 ./sample-emu.sh
	 This will run test_Emu++ with the statefile as the first argument, this creates an emulator object which is currently sampled by entering locations in the parameter space through std in. To test this type in a few locations, pressing return after each one. The code prints out the location it read and the mean and variance on new lines. When you're finished stop the program with EOF: ctrl-d or just quit with ctrl-c

3.5) use the pre-written sample locations:
		 ./sample-emu.sh < sample_locations.dat
		 	 




