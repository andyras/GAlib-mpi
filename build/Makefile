# name of main executable
### This line should be just about all you have to change ###
BIN = dynamix
#############################################################

#DEBUG_MAKE = "foo"

# directories for headers, objects and source files
INCDIR = ./include
OBJDIR = ./obj
SRCDIR = ./src
BINDIR = ../bin

# compiler
CPP = icpc
CPP = mpic++

# flags for compiler
# optimized CPPFLAGS
# Intel
#CPPFLAGS = -O0 -Wall --std=c++11 -openmp
CPPFLAGS = -O3 -fast -xHOST -no-prec-div -Wall --std=c++11 -fopenmp -mkl -no-multibyte-chars

ifdef DEBUG_MAKE
  CPPFLAGS += -pg -g -debug
endif

# non-optimized CPPFLAGS
#CPPFLAGS = -Wall -g -O0 -openmp --std=c++11

# flags for linker
INCLUDES = -I$(INCDIR) -I$(MKLROOT)/include/fftw
# includes for CVODE
INCLUDES += -I$(HOME)/bin/include -I$(HOME)/bin/include/cvode -L$(HOME)/bin/lib -I/usr/include/sundials -I/usr/include/cvode -L/usr/lib
#LDFLAGS = -lsundials_cvode -lsundials_nvecserial -lgomp -lfftw3
# Intel flags for linker
LDFLAGS = -lsundials_cvode -lsundials_nvecserial -mkl
ifdef DEBUG_MAKE
  LDFLAGS += -pg -g -debug
endif

# make will look for .cpp files in $(SRCDIR)
vpath %.cpp $(SRCDIR)

# list of source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# object files have same name as .cpp files, but with .o extension
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,obj/%.o,$(SOURCES))

# build the main executable; this should be listed first
$(BIN): $(OBJECTS)
	#$(CPP) -o $@ $^ $(LDFLAGS) $(INCLUDES)

# automatic rule for building objects
$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

# all objects depend on object directory
$(OBJECTS): | $(OBJDIR)

# make object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean cleantest build install uninstall $(BINDIR)
clean:
	rm -rf $(OBJDIR)
	rm -f $(BIN)

cleantest:
	rm -f *.out
	rm -f *.plt
	rm -f *test

build:
	make clean
	make

install: $(BIN) | $(BINDIR)
	cp dynamix ../bin

uninstall:
	rm -rf $(BINDIR)/$(BIN)

$(BINDIR):
	mkdir -p ../bin

# TODO update testing section of this file
tests: dynamix libdynamix_outputs_test.o libdynamix_input_parser_test.o
	$(CPP) -o libdynamix_outputs_test \
	 libdynamix_outputs_test.o libdynamix_outputs.o $(LDFLAGS)
	$(CPP) -o libdynamix_input_parser_test \
	 libdynamix_input_parser_test.o libdynamix_input_parser.o $(LDFLAGS)

libdynamix_outputs_test.o: libdynamix_outputs_test.cpp \
        libdynamix_outputs.cpp libdynamix_outputs.h
	$(CPP) $(CPPFLAGS) -c libdynamix_outputs_test.cpp

libdynamix_input_parser_test.o: libdynamix_input_parser_test.cpp \
        libdynamix_input_parser.cpp libdynamix_input_parser.h
	$(CPP) $(CPPFLAGS) -c libdynamix_input_parser_test.cpp
