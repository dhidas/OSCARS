NVCC = nvcc

CUDAOBJS  = $(patsubst src/%.cu,lib/%.o,$(wildcard src/*.cu))

CUDACFLAGS = -DCUDA -cudart static -std=c++11 -shared --compiler-options '-fPIC'
INCLUDE = -Iinclude

all: $(CUDAOBJS)


lib/%.o : src/%.cu
	$(NVCC) $(CUDACFLAGS) $(INCLUDE) -c $< -o $@



clean:
	rm -f $(CUDAOBJS)

