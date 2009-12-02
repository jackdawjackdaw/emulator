#!/bin/bash
# make a dist
if [ -d "r-dist" ]; then 
		rm -rf dist
		echo "removing dist dir"
fi
mkdir r-dist
echo "making new dist dir"

# copy all the files
cp -vL r-fns/*.so r-dist
cp -v r-fns/demoPlot.R r-dist
cp -v r-fns/EmuRbind.R r-dist
cp -v r-fns/testRbind.R r-dist
cp -v r-fns/plotCovReals.R r-dist
cp ./readme-dist.txt r-dist
#done

tar cvf r-dist.tar r-dist/*
gzip r-dist.tar

echo "done"
