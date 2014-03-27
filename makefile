include makevars

SRCDIR = build/src
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,build/obj/%.o,$(SOURCES))
LDFLAGS = -lsundials_cvode -lsundials_nvecserial -mkl
INCLUDES = -I./build/include -I$(MKLROOT)/include/fftw

all: lib example

lib:
	cd ga-mpi; $(MAKE)

example: example.c
	$(CXX) -g -Wall -I. -c example.c -o example.o
	$(CXX) example.o -o example -L./ga-mpi -lga-mpi

1D: 1DArray_example.cpp
	#$(CXX) $(CXXFLAGS) -g -Wall -I/opt/gcc-4.8.0/include/c++/4.8.0/ -I . -c 1DArray_example.cpp -o 1DArray_example.o -L/opt/gcc-4.8.0/lib64
	$(CXX) $(CXXFLAGS) -g -Wall -I . -c 1DArray_example.cpp -o 1DArray_example.o $(INCLUDES)
	$(CXX) 1DArray_example.o $(OBJECTS) -o 1DArray_example $(INCLUDES) -L./ga-mpi $(LDFLAGS) -lga-mpi

build/obj/dynamix.o:
	$(MAKE) -C build clean
	$(MAKE) -C build -j 24

.PHONY: clean
	rm -f ga-mpi/*.o
	rm -f ga-mpi/*.a
