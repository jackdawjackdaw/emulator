#!/bin/bash

IMAGEPATH=~/Projects/EmuPresent/images

# smash the ps directory
# if [ -d "ps" ]; then
# 		rm -rf ps
# 		echo "removing ps dir"
# fi
# echo "making new ps dir"
# mkdir ps

echo "starting R"
## run the R to make new plots
R --vanilla <presentPlot.R > /dev/null;

cd ps

## adjusted presentPlot.R to output only blank named files sp
## you can covert them to pdf without trying to trim the 
## extension off
for FILE in `ls *.ps`
do
		OUTPUT=(basename $FILE .ps)
		echo $OUTPUT
		ps2pdf $FILE $OUTPUT.pdf
done

cp *.pdf $IMAGEPATH
cp *.png $IMAGEPATH


cd ..

