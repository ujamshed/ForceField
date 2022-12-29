#
# This Makefile generates a forcefield to optimize the 3D structure of a molecule
# for MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  U Jamshed, Student
# Date Created: Dec 15, 2022

# This is the main build makefile
# This makefile can call other makefiles to build object files, libraries or executables
# at this directory level, there's nothing to build. 
# This is the root or home directory

SRC_DIR = src/
TEST_DIR = tests/
EX_DIR = examples/

all:
	cd $(SRC_DIR); make all
	cd $(TEST_DIR); make all
	cd $(EX_DIR); make all

cleanall:
	cd $(SRC_DIR); make cleanall
	cd $(TEST_DIR); make cleanall
	cd $(EX_DIR); make cleanall