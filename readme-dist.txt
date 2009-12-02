C.Coleman-Smith
cec24@phy.duke.edu
1/12/09

This is a readme for the r only "distribution" of the emulator project.

The files included are:

libEmu.so - the library which handles all of the emulation / estimation of params

libRbind.so - the library which you should load directly into R

EmuRbind.R - this contains all of the R bindings to the C code in libRbind and libEmu, the most useful  functions are callCode and setEmulatorOptions.

demoPlot.R - currently produces a side by side plot of a simple model with draws made from the emulator's distribution

plotCovReals.R - provides functions to make the draws from the emulators cov fn

testRbind.R - contains a lot of test functions which will plot the test model in different ways.

-------------------------------------------------------------------

To get started, make sure that your LD_LIBRARY_PATH is set correctly. You can check using "ldd libRbind.so" to see if it finds libEmu.so or not, you should see a few resolved libraries and "libEmu.so => not found". 

To set the correct path run (in bash)

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`

this will set your LD_LIBRARY_PATH to include the r-dist directory.


I suggest browsing through EmuRbind.R first and then running demoPlot.R and some of the test functions in testRbind.R. 
--------------------------------------------------------------------

As mentioned above the function setEmulatorOptions(covFn, alpha) is important, you need to call this *before* you run the actual emulator (using callCode forinstance). If you don't the rest of the external functions will most likely hang or produce underinfed output.

setEmulatorOptions selects which covariance function to use in the estimation/emulation process, if you select the gaussian then the alpha parameter changes the smoothness of the cov fn. Otherwise alpha does nothing. 

As listed in the source you can select from 0:gauss, 1:full matern, 2:3/2matern, 3:5/2 matern. 

You can simplify things by running setDefaultOps() instead of setEmulatorOptions(a,b) this takes no arguments and selects a gaussian with alpha=1.9

