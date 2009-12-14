#!/bin/bash

rm ps/*.pdf

## run the R to make new plots
R --vanilla <presentPlot.R> dump.txt

cd ps
for FILE in `ls *.ps`
do
		ps2pdf $FILE $FILE.pdf
done

cp *.pdf ~/Projects/EmuPresent/images

cd ..