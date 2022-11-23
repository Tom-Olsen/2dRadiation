# Compiler:
CC=g++
# Compiler Flags:
DEBUGFLAGS=#-fsanitize=address -static-libasan -g
CXXFLAGS=-std=c++17 -fopenmp #-Wfatal-errors#-ftemplate-depth=2000
POSTFLAFS=-ljsoncpp
OPTIMFLAGS=-O3 -ffast-math

SRC_DIR  := src
OBJ_DIR  := obj
SRC_FILES  := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES  := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))
INC_FILES  := $(wildcard $(SRC_DIR)/*.hh)


all: test paper

test: $(OBJ_FILES) $(INC_FILES) test.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)

paper: $(OBJ_FILES) $(INC_FILES) paper.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -c -o $@ $< $(POSTFLAFS)





#.Phony: run
#run:
#	./main 0 100 100 8
#	./main 1 100 100 8

.Phony: clean
clean:
	-rm obj/*.o
	-rm test
	-rm paper

.Phony: outputClean
outputClean:
	rm -r output/*