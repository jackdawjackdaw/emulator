#!/bin/bash

for i in {1..30} 
do
./emulator -p 2 < ./input/2d-gauss-emu-sample-10.txt; ./coverage < ./input/2d-gauss-emu-sample-10.txt;
done

for i in {1..30} 
do
./emulator -p 2 < ./input/2d-gauss-emu-sample-20.txt; ./coverage < ./input/2d-gauss-emu-sample-20.txt;
done

for i in {1..30} 
do
./emulator -p 2 < ./input/2d-gauss-emu-sample-30.txt; ./coverage < ./input/2d-gauss-emu-sample-30.txt;
done

for i in {1..30} 
do
./emulator -p 2 < ./input/2d-gauss-emu-sample-40.txt; ./coverage < ./input/2d-gauss-emu-sample-40.txt;
done

for i in {1..30} 
do
./emulator -p 2 < ./input/2d-gauss-emu-sample-50.txt; ./coverage < ./input/2d-gauss-emu-sample-50.txt;
done


./emulator -p 2 < ./input/2d-gauss-emu-sample-100.txt; ./coverage < ./input/2d-gauss-emu-sample-100.txt;
./emulator -p 2 < ./input/2d-gauss-emu-sample-100.txt; ./coverage < ./input/2d-gauss-emu-sample-100.txt;
./emulator -p 2 < ./input/2d-gauss-emu-sample-100.txt; ./coverage < ./input/2d-gauss-emu-sample-100.txt;
./emulator -p 2 < ./input/2d-gauss-emu-sample-100.txt; ./coverage < ./input/2d-gauss-emu-sample-100.txt;
./emulator -p 2 < ./input/2d-gauss-emu-sample-100.txt; ./coverage < ./input/2d-gauss-emu-sample-100.txt;