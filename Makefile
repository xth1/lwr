# Copyright: Thiago Da Silva Arruda <Thiago.xth1@gmail.com>
# All rights reserved

# LWR source files.
LWR_SRC =  src/GRIMMInterface.cpp src/Permutation.cpp src/Moviment.cpp src/Solver.cpp\
 src/DataStruct.cpp src/Util.cpp  src/InstanceGroup.cpp\
 src/SequentialDataAnalyser.cpp src/ResultsBuilder.cpp\
 src/Solvers/ExactUnitary.cpp src/Solvers/HeuristicChooseReversal.cpp\
 src/Solvers/ImproveReversals.cpp src/Solvers/GRASPHeuristic.cpp\
 src/Solvers/NonDeterministicSolution.cpp\
 src/misc/Parameter.cpp src/misc/RandomGenerator.cpp

LWR_OBJ = $(LWR_SRC:.cpp=.o)

LIB_DIR = ./lib/

#LWR
LWR_INC = -I./include/Solvers -I./include/misc -I./include/
LWR_LIB = $(LIB_DIR)liblwr.a
LWR_LIB_FLAG = -L$(LIB_DIR) -l lwr

#GRIMM
GRIMM_LIB  = $(LIB_DIR)/libgrimm.a
GRIMM_LIB_FLAG = -L$(LIB_DIR) -lgrimm
GRIMM_INC = -I./include/GRIMM/

# include directories
INCLUDES = $(GRIMM_INC) $(LWR_INC)

# library paths and flags
LIBS_FLAGS = -L$(LIB_DIR) -llwr -lgrimm -L../ -L/usr/local/lib -L/usr/lib

# C++ compiler flags 
CCFLAGS = -g -O3  -std=c++11 -Wall -Wextra -Winit-self -Wmain -Wsequence-point  -Wno-unused-result

# For some platforms could be necessery to change the compiler flag from 
#-std=c++11 to one of the following, Good luck doing That!! ;D
#{c++11,c++0x,gnu++11,gnu++0x}

# compiler
CC = g++

# Build dir
BUILD_DIR = ./build

BIN_DIR = /usr/local

#TOOLS
LWR_METAHEURISTIC_TOOL = lwr 

all: create_dirs $(GRIMM_LIB) $(LWR_SRC) $(LWR_LIB) $(LWR_METAHEURISTIC_TOOL)
.cpp.o:
	$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)

$(LWR_LIB): $(LWR_OBJ) 
	ar rcs $(LWR_LIB) $(LWR_OBJ)

$(LWR_METAHEURISTIC_TOOL):
	$(CC)  $(CCFLAGS) $(LIBS)  -o  $(BUILD_DIR)/$(LWR_METAHEURISTIC_TOOL) ./src/tools/TestMetaheuristic.cpp $(INCLUDES) $(LIBS_FLAGS)

#External libary
$(GRIMM_LIB): clean_grimm
	cd ./src/GRIMM; make; cd ../../

depend: dep

dep:
	makedepend -- $(CFLAGS) -- $(INCLUDES) $(LWR_SRC)

install:
	@mv $(BUILD_DIR)/$(LWR_METAHEURISTIC_TOOL) $(BIN_DIR)

uninstall:
	rm -rf  $(BIN_DIR)/$(LWR_METAHEURISTIC_TOOL)

# create dirs
create_dirs:
	mkdir -p ./lib; mkdir -p build

# create dirs
remove_dirs:
	rmdir ./lib; rmdir build

# Clean the 
clean: clean_grimm clean_build remove_dirs
	rm -rf $(LIB_DIR)*; rm -f $(LWR_OBJ) $(OUT)

clean_grimm:
	rm -rf $(LIB_GRIMM);cd ./src/GRIMM; make clean

clean_build:
	rm -rf $(BUILD_DIR)/*
