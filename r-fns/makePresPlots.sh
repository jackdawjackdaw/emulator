#!/bin/bash

cd ps
for FILE in `ls *.ps`
do
		ps2pdf $FILE $FILE.pdf
done

cp *.pdf ~/Projects/EmuPresent/images

cd ..