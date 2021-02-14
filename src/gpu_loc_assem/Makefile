NVFLAGS=--std=c++11
ifeq ($(DEBUG_CPU),TRUE)
	NVFLAGS+=-DDEBUG_CPU
endif

ifeq ($(DEBUG_GPU),TRUE)
	NVFLAGS+=-DDEBUG_GPU
endif

ht_loc: main.cpp driver.hpp
	nvcc -x cu $(NVFLAGS) -I. --library-path=./build --library=GPU_LOC_STATIC main.cpp -o ./ht_loc 
clean:
	rm ./ht_loc
