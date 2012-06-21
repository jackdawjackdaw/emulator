--==-- --==-- --==-- --==-- --==-- --==-- 

Tests/Examples for interactive_emulator

21.06.2012, ccs, cec24@phy.duke.edu
substantial portions of the code by Hal.Canary

--__-- --__-- --__-- --__-- --__-- --__-- 


Brief: The program interactive_emulator provides a commandline interface to the Gaussian Process (GP) emulator (libEMU), models can be 
trained with an arbitrary number of parameters and dimensions of output.

Requires: 
 GSL - The gnu scientific library 
 CMake - Kitware build tool
 R - for plotting output

Setup: 
	1) run ./build.sh in the project root to build and install the code. That's it!

Included Examples:
  - uni-simple
	- uni-2d-param
	- multi-simple
		
uni-simple:
 A 1d output model with a single parameter, this is the same model as the R example "simple-regression.R"

uni-2d-param:
 A 1d output model with two parameters, the same model as the R example "2d-regresion.R"

multi-simple:
 A sample of chemtreeN data, 3 parameters and 5 outputs.  						

				

	
								

