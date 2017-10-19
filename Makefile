CC = g++
LD = g++

NVCC = nvcc

# The following python variables MUST be updated for YOUR system if you wish to compile the c++ test routines
# ... Most people will not want to do this.  It is intended for development
PYPATH = /Library/Frameworks/Python.framework/Versions/3.6
PYINCLUDE = $(PYPATH)/include/python3.6m
PYLIBDIR = $(PYPATH)/lib
PYLIB = python3.6m

# CUDA Objects for the nvcc compiler
CUDAOBJS  = $(patsubst src/%.cu,lib/%.o,$(wildcard src/*.cu))
CUDACFLAGS = -DCUDA -cudart static -std=c++11 -shared --compiler-options '-fPIC'

# Flags for the c++ compiler
CFLAGS = -std=c++11 -fPIC -I$(PYINCLUDE) -Wall
LDFLAGS = -L$(PYLIBDIR) -l$(PYLIB) -O3

# c++ exe file names taken from what is in directories
EXECS = $(patsubst exe/%.cc,bin/%,$(wildcard exe/*.cc))
EXEOBJS  = $(patsubst exe/%.cc,lib/%.o,$(wildcard exe/*.cc))

# c++ objects based on filenames
OBJS = $(patsubst src/%.cc,lib/%.o,$(wildcard src/*.cc))

# Include the local include directory
INCLUDE = -Iinclude

# For normal compilation only compile the CUDA libraries
all: $(CUDAOBJS)


# Make "test" to compile the test execs into the bin directory
test: $(OBJS) $(EXEOBJS) $(EXECS)


lib/%.o : src/%.cu
	$(NVCC) $(CUDACFLAGS) $(INCLUDE) -c $< -o $@


lib/%.o : src/%.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : exe/%.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

bin/% : lib/%.o $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS)  $< -o $@



clean:
	rm -f $(CUDAOBJS) $(EXEOBJS) $(OBJS) $(EXECS)

