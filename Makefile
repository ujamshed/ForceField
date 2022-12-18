#
# This Makefile generates 3D integration program
# for MSSE Chem 279 - Fall Semester
# University of California, Berkeley
# 
# Creator:  U Jamshed, Student
# Date Created: Sept 22, 2022

# This is the main build makefile
# This makefile can call other makefiles to build object files, libraries or executables
# at this directory level, there's nothing to build. 
# This is the root or home directory

SRC_DIR = src/
TST_DIR = tests/

all:
	cd $(SRC_DIR); make all
	cd $(TST_DIR); make all

cleanall:
	cd $(SRC_DIR); make cleanall
	cd $(TST_DIR); make cleanall