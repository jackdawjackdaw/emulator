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
R --vanilla <presentPlot.R> dump.txt

cd ps
for FILE in `ls *.ps`
do
		ps2pdf $FILE $FILE.pdf
done

cp *.pdf ~/Projects/EmuPresent/images
cd ..

