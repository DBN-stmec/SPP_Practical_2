CC=gcc
CXX=mpic++
CFLAGS=-O3 -std=c99
CXXFLAGS=-O3
EXECS=heated-plate heated-plate-non-blocking

all: $(EXECS)

heated-plate: heated-plate-parallel.cpp
	$(CXX) -o $@ $<  $(CXXFLAGS)
heated-plate-non-blocking: heated-plate-parallel-non-blocking.cpp
	$(CXX) -o $@ $<  $(CXXFLAGS)
clean:
	rm -f $(EXECS) *.o
