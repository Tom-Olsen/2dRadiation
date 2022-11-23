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


all: main test paper eigen

main: $(OBJ_FILES) $(INC_FILES) main.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)

test: $(OBJ_FILES) $(INC_FILES) test.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)

paper: $(OBJ_FILES) $(INC_FILES) paper.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -c -o $@ $< $(POSTFLAFS)



eigen: obj/Utility.o eigen.cpp
	$(CC) $(CXXFLAGS) $(DEBUGFLAGS) $(OPTIMFLAGS) -o $@ $^ $(POSTFLAFS)





.Phony: run
run:
	./main 0 100 100 8
	./main 1 100 100 8

.Phony: clean
clean:
	-rm obj/*.o
	-rm main
	-rm test
	-rm paper
	-rm eigen

.Phony: outputClean
outputClean:
	rm -r output/*