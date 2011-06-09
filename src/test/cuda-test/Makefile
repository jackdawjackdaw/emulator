INCLUDE = -I/xtmp/cec24/NVIDIA_gpu_computing_sdk/shared/inc
LIBS = -lgsl -lgslcblas -lcublas -lcuda

approxInv: approxInvCuda.cpp
	nvcc -c   approxInvCuda.cpp 

test: approxInv
	nvcc approxInvCuda.o ../approxInverse.c -DNPOINTS=512 $(LIBS)

test2: 
	nvcc -x=c++ -o approxCuda approxInvCuda.cpp ../approxInverse.c $(INCLUDE) -DNPOINTS=512 $(LIBS)