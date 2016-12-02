CC = g++
LD = g++
NVCC = nvcc

PYVERSION = 2.7
PYPATH = /System/Library/Frameworks/Python.framework/Versions/$(PYVERSION)

CFLAGS = -DCUDA -Wall -pedantic -O3 -pthread -std=c++14 -fPIC -Wno-write-strings
CUDACFLAGS = -DCUDA
LIBS = -Llib -L$(PYPATH)/lib/python$(PYVERSION) -lpython -stdlib=libstdc++ -framework Foundation -lc++ -L/usr/local/cuda/lib -lcuda -lcudart_static
INCLUDE = -Iinclude -I$(PYPATH)/include/python$(PYVERSION)

OBJS  = $(patsubst src/%.cc,lib/%.o,$(wildcard src/*.cc))
CUDAOBJS  = $(patsubst src/%.cu,lib/%.o,$(wildcard src/*.cu))
EXECS = $(patsubst exe/%.cc,bin/%,$(wildcard exe/*.cc))
EXEOBJS  = $(patsubst exe/%.cc,lib/%.o,$(wildcard exe/*.cc))

SOLIB = lib/sr.so




all: $(OBJS) $(CUDAOBJS) $(EXEOBJS) $(EXECS) $(SOLIB)


lib/sr.so : $(OBJS) $(CUDAOBJS)
	$(LD) -shared $(LIBS) $(OBJS) $(CUDAOBJS) -o $@

lib/%.o : src/%.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : src/%.cu
	$(NVCC) $(CUDACFLAGS) $(INCLUDE) -c $< -o $@


lib/%.o : exe/%.cc
	$(CC) -Wall $(CFLAGS) $(INCLUDE) -c $< -o $@


bin/% : $(OBJS) $(CUDAOBJS) lib/%.o
	$(LD) $(LIBS) $(OBJS) $(CUDAOBJS) lib/$*.o -o bin/$*





clean:
	rm -f $(EXECS) lib/*.o $(SOLIB)

