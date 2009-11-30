#!/bin/bash
# make a dist
if [ -d "dist" ]; then 
		rm -rf dist
		echo "removing dist dir"
fi
mkdir dist
echo "making new dist dir"

# copy all the files
cp -vL r-fns/*.so dist
cp -v r-fns/demoPlot.R dist
cp -v r-fns/EmuRbind.R dist
cp -v r-fns/testRbind.R dist
cp -v r-fns/plotCovReals.R dist

#done
echo "done"
