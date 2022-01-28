CC = nvcc
PROGRAM = CUDAMCML
#CFLAGS = -arch=sm_13 --ptxas-options=-v -keep
#CFLAGS = -arch=sm_35 --ptxas-options=-v -use_fast_math #-
#CFLAGS = --ptxas-options=-v -use_fast_math #-maxrregcount=32 #-ftz=true -prec-div=false -prec-sqrt=false -maxrregcount=32
#CFLAGS = -gencode arch=compute_35,code=[sm_35,compute_35] #--ptxas-options=-v -use_fast_math

IFLAGS = -I /usr/local/cuda/include -I $(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc
LFLAGS = -L /usr/local/cuda/lib64 -L $(HOME)/NVIDIA_GPU_Computing_SDK/C/common/lib/linux -L $(HOME)/NVIDIA_GPU_Computing_SDK/C/lib -lcuda -lcudart #-lcutil
SOURCE = ./src/CUDAMCMLmain.cu
DEP = ./src/CUDAMCMLmain.cu ./src/CUDAMCMLio.cu ./src/CUDAMCMLmem.cu ./src/CUDAMCMLrng.cu ./src/CUDAMCMLtransport.cu ./src/CUDAMCML.h

$(PROGRAM): $(DEP)
	$(CC) $(CFLAGS) $(SOURCE) $(IFLAGS) $(LFLAGS) -o $(PROGRAM)

clean:
	\rm -f *~ *.o *.linkinfo CUDAMCML



