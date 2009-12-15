#!/bin/bash

# smash the ps directory
if [ -d "ps" ]; then
		rm -rf ps
		echo "removing ps dir"
fi
echo "making new ps dir"
mkdir ps

echo "starting R"
## run the R to make new plots
R --vanilla <presentPlot.R > /dev/null;

cd ps

## adjusted presentPlot.R to output only blank named files sp
## you can covert them to pdf without trying to trim the 
## extension off
for FILE in `ls *.`
do
		ps2pdf $FILE $FILE.pdf
done

cp *.pdf ~/Projects/EmuPresent/images

cd ..

