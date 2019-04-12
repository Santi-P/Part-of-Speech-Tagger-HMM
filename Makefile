CXX		  := g++
CXX_FLAGS :=  -O2 -std=c++17

BIN		:= bin
SRC		:= src
INCLUDE	:= include
LIB		:= lib
DATA	:= data

LIBRARIES	:=
EXECUTABLE	:= hmm


all: $(BIN)/$(EXECUTABLE)

install: 
	

tests:
	#$(BIN)/$(EXECUTABLE)  data/tiger_spaces.txt data/tiger_spaces.txt
	$(BIN)/$(EXECUTABLE)  data/tiger_train.txt data/tiger_test.txt > results.txt
	python -c 'from eval import evaluate; evaluate("data/tiger_test.txt", "results.txt")'




documentation:Doxyfile
	doxygen



run:
	
	clean all
	clear
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)
	

clean:
	-rm $(BIN)/*
	-rm -r doc/*
	-rm results.txt
