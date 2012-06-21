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

Interactive_emulator command line args:
    interactive_emulator estimate_thetas INPUT_MODEL_FILE MODEL_SNAPSHOT_FILE [OPTIONS]
  or
    interactive_emulator interactive_mode MODEL_SNAPSHOT_FILE

  Most importantly the first argument should be either: "estimate_thetas" or "interactive_mode"
	estimate_thetas will train an emulator on the given data set producing a state file
	interactive_mode reads parameter space locations from stdin and outputs the emulator mean + variance + covariance matrix (coming soon) at these points.

	Options for estimate_thetas can include:
    --regression_order=0
    --regression_order=1
    --regression_order=2
    --regression_order=3
    --covariance_fn=POWER_EXPONENTIAL
    --covariance_fn=MATERN32
    --covariance_fn=MATERN52
  The defaults are regression_order=0 and covariance_fn=POWER_EXPONENTIAL.


Interactive_emulator input data Format:
		 Input files should be laid out as follows

INPUT_MODEL_FILE (With multi-output) FORMAT:
  Models with multivariate output values y = y_1...y_{number_outputs} which we will think of as rows of a matrix, in the following spec
  number_outputs = t
  BEGIN EXAMPLE
    number_outputs
    number_params 
    number_model_points
    X[0,0]
    ...
    X[0,number_params-1]
    X[1,0]
    ...
    X[number_model_points-1,number_params-1]
    Y[0, 0]
    Y[0, 1]
    ...
    Y[0, number_outputs-1]
    ...
    Y[1, 1]
    ... 
    Y[number_model_points-1, number_outputs-1]
   END EXAMPLE

  number_outputs, number_params and number_model_points should be
  positive integers.  X[i,j] and Y[i,j] will be read as
  double-precision floats.

  You don't have to use newlines to separate values.  If fscanf(fp,
  "%lf%*c", ptr) can read the value, it will work.
	
				

	
								

