include makevars

all: lib example

lib:
	cd ga-mpi; $(MAKE)

example: example.c
	$(CXX) -g -Wall -I. -c example.c -o example.o
	$(CXX) example.o -o example -L./ga-mpi -lga-mpi

1D: 1DArray_example.cpp
	$(CXX) $(CXXFLAGS) -g -Wall -I/opt/gcc-4.8.0/include/c++/4.8.0/ -I . -c 1DArray_example.cpp -o 1DArray_example.o -L/opt/gcc-4.8.0/lib64
	$(CXX) 1DArray_example.o -o 1DArray_example -L./ga-mpi -lga-mpi

.PHONY: clean
	rm -f ga-mpi/*.o
	rm -f ga-mpi/*.a
