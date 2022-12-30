#
# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  U Jamshed, Student
# Date Created: Dec 15nd, 2022
#
# This makefile creates object files containing tests that are executable (i.e. they have main() functions)
# They test the functionalities in ROOT_DIRECTORY. After creating the executable files they
# are relocated to the ../bin directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -O3 -o 

# The directory where the include files needed to create the library objects are
INC_DIR = ../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../lib
BIN_DIR     = ../bin

# The next line contains the list of tests created by this Makefile.
EXECUTABLES = unit_tests

unit_tests: $(INC_DIR)/utility.h $(INC_DIR)/input.h $(INC_DIR)/energy.h
	$(CXX) $(CXXFLAGS) unit_tests unit_tests.cpp -I$(INC_DIR) $(LIB_DIR)/utility.o $(LIB_DIR)/input.o $(LIB_DIR)/energy.o $(LIB_DIR)/optimizers.o $(LIB_DIR)/output.o -larmadillo

all: $(EXECUTABLES)

cleanall:
	rm $(EXECUTABLES)