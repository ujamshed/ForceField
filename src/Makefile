#
# MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  U Jamshed, Student
# Date Created: Sept 22nd, 2022
#
# This makefile creates objects from the src directory of functions and procedures in C++
# after creation the objects are relocated to the ../lib directory

# GNU C++ Compiler
CXX         = g++      

# compiler flags -g debug, -O3 optimized version -c create a library object
CXXFLAGS    = -O3 -c      

# The directory where the include files needed to create the library objects are
INC_DIR = ../include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../lib

# The next line contains the list of math object files created by this Makefile.
# they are the targets this makefile builds
DATA_STRUCT_OBJS = utility.o input.o energy.o 

utility.o: $(INC_DIR)/utility.h
	$(CXX) $(CXXFLAGS) utility.cpp -I$(INC_DIR)
	mv utility.o $(LIB_DIR)

input.o: $(INC_DIR)/input.h
	$(CXX) $(CXXFLAGS) input.cpp -I$(INC_DIR)
	mv input.o $(LIB_DIR)

energy.o: $(INC_DIR)/energy.h
	$(CXX) $(CXXFLAGS) energy.cpp -I$(INC_DIR)
	mv energy.o $(LIB_DIR)

all: $(DATA_STRUCT_OBJS)

cleanall:
	cd $(LIB_DIR); rm $(DATA_STRUCT_OBJS)