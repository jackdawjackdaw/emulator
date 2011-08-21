#!/bin/zsh

sizes=( 32 64 128 256 512 1024 2048 4096 8192 16384 32768)


for i in $sizes; do
		echo "running $i";
		#gcc approxInverse.c -o invApprox -lgsl -lgslcblas -lm -O3 -g -DNPOINTS=$i;
		nvcc -x=c++ -o approxCuda approxInvCuda.cpp ../approxInverse.c -DNPOINTS=$i -I/xtmp/cec24/NVIDIA_gpu_computing_sdk/shared/inc  -lgsl -lgslcblas -lcublas -lcuda > /dev/null

		./approxCuda;
done
