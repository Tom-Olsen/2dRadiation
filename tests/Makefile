# Compiler:
CC=g++
# Compiler Flags:
DEBUGFLAGS=#-fsanitize=address -static-libasan -g
CXXFLAGS=-std=c++17 -fopenmp #-Wfatal-errors#-ftemplate-depth=2000
POSTFLAFS=-ljsoncpp
OPTIMFLAGS=-O3 #-ffast-math

SRC_DIR  := ../src
OBJ_DIR  := ../obj
SRC_FILES  := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES  := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))
INC_FILES  := $(wildcard $(SRC_DIR)/*.hh)



all: testStencil #testMetric2D

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -c -o $@ $< $(POSTFLAFS)

testStencil: $(OBJ_FILES) $(INC_FILES) testStencil.cpp 
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)
	
testMetric2D: $(OBJ_FILES) $(INC_FILES) testMetric2D.cpp 
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)



.Phony: clean
clean:
	-rm testStencil
	-rm testMetric2D