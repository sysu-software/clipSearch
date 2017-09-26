
# ==========================
# clipSearch Makefile
# (c) 2010 Jian-Hua Yang
# yangjh7@mail.sysu.edu.cn
# ==========================
# define our object and binary directories
BIN_DIR	 = bin
CXX		 = g++
CC       = gcc
CXXFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -fPIC $(INCLUDE)
LIBS	 = -lz -lm
RNA_DIR  = thirdUtils/RNAfoldLib
UTILS_DIR = bioUtils
PROGRAM_DIR = src

all:
	cd $(RNA_DIR); make
	cd $(UTILS_DIR); make
	cd $(PROGRAM_DIR)/clipSearch; make

clean:
	cd $(RNA_DIR); make clean
	cd $(UTILS_DIR); make clean
	cd $(PROGRAM_DIR)/clipSearch; make clean
	cd $(BIN_DIR); rm clipSearch
	
test:
	@cd test_data; bash run_test.sh
	