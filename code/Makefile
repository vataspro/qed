CXX=g++
CXXFLAGS=-std=c++17

SRCS = main.cpp lattice.cpp ranlxs.cpp mc.cpp

all:: ranlxs.o lattice.o mc.o
	$(CXX) $(CXXFLAGS) $(SRCS) -lm -o exec

clean:
	rm *.o exec
